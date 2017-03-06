-- |
-- Module      : Math.FunctionalAnalysis.L2Function.R1
-- Copyright   : (c) Justus Sagemüller 2017
-- License     : GPL v3
-- 
-- Maintainer  : (@) sagemueller $ geo.uni-koeln.de
-- Stability   : experimental
-- Portability : portable
-- 
{-# LANGUAGE TypeFamilies               #-}
{-# LANGUAGE ConstraintKinds            #-}
{-# LANGUAGE ScopedTypeVariables        #-}
{-# LANGUAGE UnicodeSyntax              #-}
{-# LANGUAGE FlexibleContexts           #-}

module Math.FunctionalAnalysis.L2Function.R1 (
                         UnitL2(..)
                       , toUniformSampled
                       , fromUniformSampled
                       , SampleMode(..)
                       , SigSampleConfig(..)
            ) where

import qualified Data.Map as Map
import qualified Data.Vector.Storable as UArr
import qualified Data.Vector as BArr
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as Arr
import qualified Numeric.FFT.Vector.Plan as FFT
import qualified Numeric.FFT.Vector.Invertible as FFT

import Control.Arrow
import Data.Monoid ((<>))
import Data.Complex

import Data.Foldable (fold)

nSubDivs :: Int
nSubDivs = 8

subdivFreq :: Double
subdivFreq = 1 / fromIntegral nSubDivs

data SampleMode = OutlappingDFT

type UABSample c = (Integral c, Bounded c, UArr.Storable c)

-- | A square-integrable, bandlimited function with compact support
--   on the unit interval. 
data UnitL2 c s = UnitL2 {
    unitL2AmplitudeFactor :: !s
  , unitL2ExtremeSampleVals :: !(c,c)
  , unitL2LoFreqSamples :: UArr.Vector (Complex c)
  , unitL2LoFreqSampleMode :: !SampleMode
  , unitL2Subdivisions :: BArr.Vector (UnitL2 c s)
  }

data HomogenSampled v = HomogenSampled {
       _homogenSampledRange :: (HIndex v, HIndex v)
     , _values :: UArr.Vector v      -- ^ Starts and ends with a zero, which is not
                                     --   considered as part of the covered range.
                                     --   Index position 0 is mid between this
                                     --   “ghost cell” and the first value cell.
                                     --   Each value thus represents one /interval/
                                     --   of size ℓ/𝑛, where 𝑛+2 is the length of
                                     --   the values array
     , _derivatives :: UArr.Vector v -- ^ Corresponds to function as sampled on [0,1].
                                     --   This is also zero-bounded.
     } deriving (Show)
homogenSampled :: (Fractional v, UArr.Storable v)
                    => (HIndex v, HIndex v) -> UArr.Vector v -> HomogenSampled v
homogenSampled range vs | n > 1
           = HomogenSampled range bvs derivs
 where derivs = Arr.cons 0 . (`Arr.snoc`0)
                  $ Arr.imap (\j _ -> let vp = Arr.unsafeIndex bvs j
                                          vn = Arr.unsafeIndex bvs (j+2)
                                      in (vn - vp)/(2*h) )
                             vs
       n = Arr.length vs
       h = 1 / fromIntegral n
       bvs = Arr.cons 0 $ Arr.snoc vs 0

unitHomogenSampled :: (Fractional v, UArr.Storable v)
                      => UArr.Vector v -> HomogenSampled v
unitHomogenSampled = homogenSampled (0,1)

homogenSampledInformation :: (RealFrac v, UArr.Storable v) => HomogenSampled v -> Int
homogenSampledInformation (HomogenSampled (start,end) vs _) = ceiling $ l * ν
 where ν = fromIntegral $ Arr.length vs - 2
       l = end - start

class DynamicDimension a where
  dynDimension :: a -> Int
instance UArr.Storable v => DynamicDimension (UArr.Vector v) where
  dynDimension = Arr.length
instance (RealFrac v, UArr.Storable v) => DynamicDimension (HomogenSampled v) where
  dynDimension = homogenSampledInformation

subdivideHomogenSampled :: Fractional v 
              => Int -> HomogenSampled v -> BArr.Vector (HomogenSampled v)
subdivideHomogenSampled n (HomogenSampled (start,end) vs dvs)
       = Arr.generate n (\j -> let rStart = start + fromIntegral j * lr
                                   rEnd = rStart + lr
                               in HomogenSampled (rStart,rEnd) vs dvs )
 where l = end - start
       lr = l / fromIntegral n


-- | 0 corresponds to the first nonzero entry of a 'HomogenSampled', 1 to the last.
type HIndex v = v

resampleHomogen :: (RealFrac v, UArr.Storable v, Show v)
            => (HIndex v, HIndex v) -> Int -> HomogenSampled v -> HomogenSampled v
resampleHomogen (start', end') n (HomogenSampled (start₀, end₀) vs dvs)
     | end > start, n > 0, nOld > 0
          = HomogenSampled (0,1)
                           (preZeroes<>vals<>postZeroes)
                           (preZeroes<>derivs<>postZeroes)
     | otherwise  = HomogenSampled (0,1) (Arr.replicate 2 0) (Arr.replicate 2 0)
 where -- see
       -- https://raw.githubusercontent.com/leftaroundabout/sobolev-spaces/master/derivation/interpolation-alignment.svg
       [nPreZeroes, nPostZeroes] = max 0 . min (n+1) . round
                              <$> [(-start-hOld/2) / h, (end-1-hOld/2) / h]
       [preZeroes, postZeroes] = (`Arr.replicate`0) . (+1) <$> [nPreZeroes, nPostZeroes]
       nNonzero = n - nPreZeroes - nPostZeroes
       h = l / fromIntegral n
       ts = Arr.generate nNonzero $
              \j -> start + h * (fromIntegral (nPreZeroes + j) + 1/2)
       vals = Arr.map (\t -> let jp = max 0 . min nOld . round $ t/hOld
                                 y₀ = vs Arr.! jp; y₁ = vs Arr.! (jp+1)
                                 ð₀ = (dvs Arr.! jp)*hOld; ð₁ = (dvs Arr.! (jp+1))*hOld
                                 ξ = t / hOld - fromIntegral (jp-1) - 1/2
                             in y₀
                                + ξ * (ð₀
                                      + ξ * (3*(y₁ - y₀) - 2*ð₀ - ð₁
                                            + ξ * (ð₁ + ð₀ + 2*(y₀ - y₁))))
                      ) ts
       derivs = Arr.map (\t -> let jp = max 0 . min nOld . round $ t/hOld
                                   ð₀ = (dvs Arr.! jp)*hOld; ð₁ = (dvs Arr.! (jp+1))*hOld
                                   ξ = t / hOld - fromIntegral (jp-1)
                               in (ð₀ + ξ*(ð₁-ð₀)) / h
                        ) ts
       nOld = Arr.length vs - 2
       hOld = 1 / fromIntegral nOld
       l = end - start
       start = start₀ + l₀ * start'
       end = start₀ + l₀ * end'
       l₀ = end₀ - start₀

{-# SPECIALISE resampleHomogen :: (Double, Double) -> Int -> HomogenSampled Double -> HomogenSampled Double #-}

      
-- | Efficient only when partially applied to data array
--   and then evaluated at multiple x-points 
lfCubicSplineEval :: UArr.Vector Double -> Double -> Double
lfCubicSplineEval ys
  | n > 1      = lcse
  | n < 1      = const 0
  | otherwise  = \x
     -> if x<0 || x>1 then 0
                      else ys!0 * ((2*x-1)^2 - 1)^2
 where lcse x | i < 0      = 0
              | i < 1      = ξ^2 * (3*y₁ - ð₁ + ξ*(ð₁ - 2*y₁))
              | i < n      = y₀
                           + ξ * (ð₀
                                  + ξ * (3*(y₁ - y₀) - 2*ð₀ - ð₁
                                         + ξ * (ð₁ + ð₀ + 2*(y₀ - y₁))))
              | i == n     = (1-ξ)^2 * (3*y₀ + ð₀ + (1-ξ)*(-ð₀ - 2*y₀))
              | otherwise  = 0
        where hx = x / δx
              i = floor hx
              ξ = hx - fromIntegral i
              y₀ = ys!(i-1)
              y₁ = ys!i
              ð₀ = 0.5 * δys!(i-1)
              ð₁ = 0.5 * δys!i
       δys = UArr.imap (\i _
                  -> if i<1
                      then ys!1
                      else if i<n-1
                            then ys!(i+1) - ys!(i-1)
                            else -ys!(n-2)
              ) ys
       n = Arr.length ys
       δx = 1 / fromIntegral (n+1)

cubicResample :: Arr.Vector v Double => Int -> v Double -> v Double
cubicResample n ys = Arr.generate n $ \i -> spline $ fromIntegral (i+1) / fromIntegral (n+1)
 where spline = lfCubicSplineEval $ Arr.convert ys

fourierTrafo :: UArr.Vector (Complex Double) -> HomogenSampled Double
fourierTrafo = prepareTrafos 2 FFT.dft process
 where process n' fft = homogenSampled (0.25, 0.75)
                      . Arr.zipWith (\μ -> realPart . (μ*)) μs
                      . FFT.execute fft
                      . \αs -> Arr.zipWith (*) ηs αs <> Arr.replicate (n' - Arr.length αs) 0
        where μs = Arr.generate n' $ \j -> let t = (1-n)/n + 2*fromIntegral j/n
                                           in cis $ -pi * t/2
              ηs = Arr.generate n' $ \k -> cis (-pi*fromIntegral k*(1-n)/n)
              n = fromIntegral n'

invFourierTrafo :: HomogenSampled Double -> UArr.Vector (Complex Double)
invFourierTrafo = prepareTrafos 2 FFT.idft process
 where process n' ift = Arr.zipWith (*) ηs
                      . FFT.execute ift
                      . (\(HomogenSampled (0,1) vs _)
                        -> Arr.zipWith (\(rμ:+iμ) z -> rμ*z :+ iμ*z) μs
                            $ Arr.slice 1 n' vs )
                      . resampleHomogen (-0.5, 1.5) n'
        where μs = Arr.generate n' $ \j -> let t = (1-n)/n + 2*fromIntegral j/n
                                           in cis $ pi * t/2
              ηs = Arr.generate (n'`div`2) $ \k -> 2 * cis (pi*fromIntegral k*(1-n)/n)
              n = fromIntegral n'

prepareTrafos :: (UArr.Storable a, UArr.Storable b, DynamicDimension c)
                    => Int -> FFT.Transform a b
                     -> (Int -> FFT.Plan a b -> c -> r)
                     -> c -> r
prepareTrafos sizeFactor transform env = \αs
  -> let lookupTrafo nmin = case Map.lookupGE nmin plans of
          Just (_,f) -> f αs
          Nothing -> lookupTrafo $ (nmin*2)`div`3
     in lookupTrafo $ dynDimension αs * sizeFactor
 where plans = Map.fromList [ (n, env n $ FFT.plan transform n)
                            | p <- [3..10]
                            , υ <- [1,0]
                            , let n = 2^p - υ*2^(p-2) -- [6,8,12,16,24,32,48.. 1024]
                            ]



data SigSampleConfig = SigSampleConfig {
      _maxFFTSize :: Int  -- ^ Not considered at the moment; max FFT is always 1024. 
    , _infoPerStage :: Int
    , _maxLocalInfo :: Int
    , _longrangeBandwidth :: Double
    , _noiseFloorLevel :: Double
    }

-- | A phase-constant second-order filter, using one first-order pass
--   in each direction.
simpleIIRLowpass, simpleIIRHighpass
                 :: Double       -- ^ Cutoff frequency, normalised to 1/τ
                                 --   (total data length, not sampling frequency)
                 -> UArr.Vector Double
                 -> UArr.Vector Double
simpleIIRLowpass ω ys = Arr.postscanr' (\y carry -> (1-η)*carry + η*y) 0
                      $ Arr.postscanl' (\carry y -> (1-η)*carry + η*y) 0 ys
 where η = ω / (ω + fromIntegral (Arr.length ys))

simpleIIRHighpass ω ys = Arr.zipWith (-) ys $ simpleIIRLowpass ω ys

toUniformSampledLike :: UABSample c
           => HomogenSampled Double -> UnitL2 c Double -> UArr.Vector Double
toUniformSampledLike range@(HomogenSampled (start,end) vs _)
                        (UnitL2 μLr _ quantisedLr _ subchunks) = result
 where -- see
       -- https://raw.githubusercontent.com/leftaroundabout/sobolev-spaces/master/derivation/sampling-alignment.svg
       nRef = Arr.length vs - 2
       ℓ = end - start
       hRef = 1 / fromIntegral nRef
       h = hRef / ℓ
       t₀Ref = - start / ℓ
       nBefore = round $ start / hRef
       iEnd = round $ end / hRef
       t₀ = t₀Ref + fromIntegral nBefore * h
       tEnd = t₀Ref + fromIntegral iEnd * h
       n = iEnd - nBefore
       backTransformed = fourierTrafo
                         $ Arr.map ((*realToFrac μLr) . fromIntegralℂ) quantisedLr
       HomogenSampled _ resultLr' _ = resampleHomogen (t₀,tEnd) n backTransformed
       result, resultLr, subResult :: UArr.Vector Double
       resultLr = Arr.slice 1 n resultLr'
       result = Arr.zipWith (+) resultLr (subResult
                             <> Arr.replicate (Arr.length resultLr
                                                    - Arr.length subResult) 0)
       subResult = fold $ Arr.zipWith toUniformSampledLike
                     (subdivideHomogenSampled nSubDivs range) subchunks

toUniformSampled :: UABSample c => Int -> UnitL2 c Double -> UArr.Vector Double
toUniformSampled n = toUniformSampledLike . unitHomogenSampled $ UArr.replicate n 0

fromUniformSampled :: ∀ c . UABSample c
        => SigSampleConfig
        -> UArr.Vector Double -> UnitL2 c Double
fromUniformSampled cfg allYs = result
 where result = chunkFromUniform cfg residuals
       residuals = residualLayers cfg allYs $ pure result

onlyLongrange :: UnitL2 c Double -> UnitL2 c Double
onlyLongrange f = f {unitL2Subdivisions = Arr.empty}

residualLayers :: UABSample c => SigSampleConfig -> UArr.Vector Double
                      -> BArr.Vector (UnitL2 c Double) -> [HomogenSampled Double]
residualLayers cfg@(SigSampleConfig _ _ _ lrBandwidth _) allYs modChunks
       = lowpassed : residualLayers cfg topResidual subchunks
 where longrange = fold $ Arr.zipWith toUniformSampledLike
              (subdivideHomogenSampled (Arr.length modChunks) lowpassed)
              (onlyLongrange <$> modChunks)
       lowpassed = unitHomogenSampled (simpleIIRLowpass lrBandwidth allYs)
       topResidual = Arr.zipWith (-) allYs longrange
       subchunks = modChunks >>= unitL2Subdivisions

chunkFromUniform :: ∀ c . UABSample c
               => SigSampleConfig
               -> [HomogenSampled Double] -- ^ Lazy list of signal parts. The head
                                          --   contains just the lowpassed version
                                          --   of the full signal we try to describe;
                                          --   the subsequent layers contain residual
                                          --   information not yet captured by that
                                          --   long-range description.
               -> UnitL2 c Double
chunkFromUniform cfg@(SigSampleConfig nChunkMax
                           infoPerStage localInfo lrBandwidth noiseLvl)
                      (lowpassed:residuals)
            = UnitL2 μLr (round $ maxNegAmplitudeLr/μLr, round $ maxPosAmplitudeLr/μLr)
                  quantisedLr OutlappingDFT
                  subchunks
 where transformedLr = invFourierTrafo lowpassed
       maxPosAmplitudeLr = Arr.foldl' maxℂ 0 transformedLr
       maxNegAmplitudeLr = Arr.foldl' minℂ 0 transformedLr
       μLr = maximum [maxPosAmplitudeLr, -maxNegAmplitudeLr, noiseLvl]
               / maxAllowedVal
       powerSpectrumLr = simpleIIRHighpass (fromIntegral nTot / subdivFreq)
                                $ Arr.map abs²ℂ transformedLr
       subdivideFurther :: Bool
       subdivideFurther = nTot > max (3*nSubDivs) localInfo
       quantisedLr :: UArr.Vector (Complex c)
       quantisedLr = Arr.map (roundℂ . (/realToFrac μLr))
            $ if subdivideFurther
               then filterFirstNPositivesOf powerSpectrumLr transformedLr
               else transformedLr
       subResiduals = transposeV
               $ map (subdivideHomogenSampled nSubDivs) residuals
       subchunks
        | subdivideFurther  = Arr.map (chunkFromUniform cfg) subResiduals
        | otherwise         = Arr.empty
       maxAllowedVal = fromIntegral (maxBound :: c) / 8
       nTot = dynDimension lowpassed
       nSingleChunk = fromIntegral nTot / subdivFreq
       nFullChunk = ceiling $ fromIntegral nTot * subdivFreq
       filterFirstNPositivesOf :: (Num n, UArr.Storable n) => UArr.Vector Double
                                   -> UArr.Vector n -> UArr.Vector n
       filterFirstNPositivesOf spect d = (`Arr.unfoldr`(0,infoPerStage))
           $ \(i,capacity) -> if capacity>0 && i<Arr.length d
                               then Just $ if spect!i > 0
                                            then (d!i, (i+1, capacity-1))
                                            else (0, (i+1, capacity))
                               else Nothing


transposeV :: (Arr.Vector v a, Arr.Vector v [a]) => [v a] -> v [a]
transposeV (h:q) = Arr.imap (\i xh -> xh : map (Arr.!i) q) h
transposeV [] = Arr.empty


fromIntegralℂ :: (Num s, Integral c) => Complex c -> Complex s
fromIntegralℂ (r:+i) = fromIntegral r :+ fromIntegral i

roundℂ :: (RealFrac s, Integral c) => Complex s -> Complex c
roundℂ (r:+i) = round r :+ round i

minℂ, maxℂ :: Ord s => s -> Complex s -> s
minℂ p (r:+i) = min p (min r i)
maxℂ p (r:+i) = max p (max r i)

abs²ℂ :: Num s => Complex s -> s
abs²ℂ (r:+i) = r^2 + i^2
