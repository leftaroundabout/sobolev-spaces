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

subdivisionsSizeFactor :: Int
subdivisionsSizeFactor = 8

nSubDivs :: Int
nSubDivs = 9

-- ℓ = 1/subdivisionsSizeFactor
-- τ = 1/subdivFreq
-- (𝑛−1) ⋅ τ + ℓ = 1
-- τ = (1 − ℓ)/(𝑛 − 1)
subdivFreq :: Double
subdivFreq = (fromIntegral nSubDivs - 1) / (1 - 1/fromIntegral subdivisionsSizeFactor)

subdivsOverlap :: Double
subdivsOverlap = 1 - fromIntegral subdivisionsSizeFactor / subdivFreq

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
       h = 1 / (fromIntegral n-1)
       bvs = Arr.cons 0 $ Arr.snoc vs 0


homogenSampledInformation :: (RealFrac v, UArr.Storable v) => HomogenSampled v -> Int
homogenSampledInformation (HomogenSampled (start,end) vs _) = ceiling $ l * ν
 where ν = fromIntegral $ Arr.length vs - 3
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

resampleHomogen :: (RealFrac v, UArr.Storable v)
            => HomogenSampled v -> (HIndex v, HIndex v) -> Int -> HomogenSampled v
resampleHomogen (HomogenSampled (start₀, end₀) vs dvs) (start', end') n
     | end > start, n > 1, nOld > 1
          = HomogenSampled (0,1)
                           (preZeroes<>vals<>postZeroes)
                           (preZeroes<>derivs<>postZeroes)
 where -- see
       -- https://raw.githubusercontent.com/leftaroundabout/sobolev-spaces/master/derivation/interpolation-alignment.svg
       [nPreZeroes, nPostZeroes] = max 0 . min (n+1) . ceiling
                              <$> [(-start-hOld) / h, (end-1-hOld) / h]
       [preZeroes, postZeroes] = (`Arr.replicate`0) . (+1) <$> [nPreZeroes, nPostZeroes]
       nNonzero = n - nPreZeroes - nPostZeroes
       h = l / (fromIntegral n - 1)
       ts = Arr.generate nNonzero $
              \j -> start + h * fromIntegral (nPreZeroes + j)
       vals = Arr.map (\t -> let jp = max 0 . min nOld . ceiling $ t/hOld
                                 y₀ = vs Arr.! jp; y₁ = vs Arr.! (jp+1)
                                 ð₀ = (dvs Arr.! jp)*hOld; ð₁ = (dvs Arr.! (jp+1))*hOld
                                 ξ = t / hOld - fromIntegral (jp-1)
                             in y₀
                                + ξ * (ð₀
                                      + ξ * (3*(y₁ - y₀) - 2*ð₀ - ð₁
                                            + ξ * (ð₁ + ð₀ + 2*(y₀ - y₁))))
                      ) ts
       derivs = Arr.map (\t -> let jp = max 0 . min nOld . ceiling $ t/hOld
                                   ð₀ = (dvs Arr.! jp)*hOld; ð₁ = (dvs Arr.! (jp+1))*hOld
                                   ξ = t / hOld - fromIntegral (jp-1)
                               in (ð₀ + ξ*(ð₁-ð₀)) / h
                        ) ts
       nOld = Arr.length vs - 2
       hOld = 1 / (fromIntegral nOld - 1)
       l = end - start
       start = start₀ + l₀ * start'
       end = start₀ + l₀ * end'
       l₀ = end₀ - start₀

{-# SPECIALISE resampleHomogen :: HomogenSampled Double -> (Double, Double) -> Int -> HomogenSampled Double #-}

      
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

fourierTrafo :: UArr.Vector (Complex Double) -> UArr.Vector Double
fourierTrafo = prepareTrafos 2 FFT.dft $ prepare &&& untwirl
 where untwirl n' = Arr.zipWith (\μ -> realPart . (μ*)) μs
        where μs = Arr.generate n' $ \j -> let t = (1-n)/n + 2*fromIntegral j/n
                                           in cis $ -pi * t/2
              n = fromIntegral n'
       prepare n' = \αs -> Arr.zipWith (*) ηs αs <> Arr.replicate (n' - Arr.length αs) 0
        where ηs = Arr.generate n' $ \k -> cis (-pi*fromIntegral k*(1-n)/n)
              n = fromIntegral n'

invFourierTrafo :: UArr.Vector Double -> UArr.Vector (Complex Double)
invFourierTrafo = prepareTrafos 1 FFT.idft
                    $ \n -> (pretwirl n . cubicResample n, postwirl n)
 where pretwirl n' = Arr.zipWith (\(rμ:+iμ) z -> rμ*z :+ iμ*z) μs
        where μs = Arr.generate n' $ \j -> let t = (1-n)/n + 2*fromIntegral j/n
                                           in cis $ pi * t/2
              n = fromIntegral n'
       postwirl n' = Arr.zipWith (*) ηs
        where ηs = Arr.generate (n'`div`2) $ \k -> 2 * cis (pi*fromIntegral k*(1-n)/n)
              n = fromIntegral n'

prepareTrafos :: (UArr.Storable a, UArr.Storable b, DynamicDimension c)
                    => Int -> FFT.Transform a b
                     -> (Int -> (c -> UArr.Vector a, UArr.Vector b -> r))
                     -> c -> r
prepareTrafos sizeFactor transform env = \αs
   -> let nmin = dynDimension αs * sizeFactor
      in case Map.lookupGE nmin plans of
          Just (n,(plan,(preproc,postproc)))
              -> postproc . FFT.execute plan $ preproc αs
 where plans = Map.fromList [ (n, (FFT.plan transform n, env n))
                            | p <- [3..10]
                            , υ <- [0,1]
                            , let n = 2^p + υ*2^(p-1) -- [8,12,16,24,32,48..]
                            ]

evalUnitL2 :: (Integral c, UArr.Storable c) => UnitL2 c Double -> Double -> Double
evalUnitL2 (UnitL2 μ _ lf OutlappingDFT subdivs) = evalAt
 where evalAt x
         | x < 0 || x > 1
                      = 0
         | nsd == 0   = lfEval x
         | i < 1 || ξ > subdivsOverlap && i < nsd
                      = lfEval x + (subdivEval!i) ξ
         | i < nsd    = lfEval x + (subdivEval!i) ξ
                                 + (subdivEval!(i-1)) (ξ + 1 - subdivsOverlap)
         | i == nsd   = lfEval x + (subdivEval!(i-1)) (ξ + 1 - subdivsOverlap)
         | otherwise  = lfEval x
        where i = floor $ x * subdivFreq
              xr = x - fromIntegral i / subdivFreq
              ξ = xr * fromIntegral subdivisionsSizeFactor
       lfEval | Arr.length lf > 0  = lfCubicSplineEval
                  (fourierTrafo $ UArr.map ((*realToFrac μ) . fromIntegralℂ) lf)
              | otherwise          = const 0
       subdivEval = Arr.map evalUnitL2 subdivs
       nsd = Arr.length subdivs


data SigSampleConfig = SigSampleConfig {
      _maxFFTSize :: Int
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
toUniformSampledLike (HomogenSampled (start,end) vs _) (UnitL2 μLr _ quantisedLr _ _)
             = Arr.slice 1 n result
 where nRef = Arr.length vs - 3
       l = end - start
       hRef = 1 / fromIntegral nRef
       h = l / fromIntegral nRef
       t₀Ref = - start / l 
       i₀ = ceiling $ start / hRef
       iEnd = ceiling $ end / hRef
       t₀ = t₀Ref + fromIntegral i₀ * h
       tEnd = t₀Ref + fromIntegral (iEnd-1) * h
       n = iEnd - i₀
       HomogenSampled _ result _ = resampleHomogen backTransformed (t₀,tEnd) n
       backTransformed = homogenSampled (0,1) . fourierTrafo
                         $ Arr.map ((*realToFrac μLr) . fromIntegralℂ) quantisedLr

toUniformSampled :: UABSample c => Int -> UnitL2 c Double -> UArr.Vector Double
toUniformSampled n = toUniformSampledLike . homogenSampled (0,1) $ UArr.replicate (n+1) 0

fromUniformSampled :: ∀ c . UABSample c
        => SigSampleConfig
        -> UArr.Vector Double -> UnitL2 c Double
fromUniformSampled cfg allYs = result
 where result = chunkFromUniform cfg residuals
       residuals = residualLayers cfg allYs $ pure result

residualLayers :: UABSample c => SigSampleConfig -> UArr.Vector Double
                      -> BArr.Vector (UnitL2 c Double) -> [HomogenSampled Double]
residualLayers cfg@(SigSampleConfig _ _ _ lrBandwidth _) allYs modChunks
       = lowpassed : residualLayers cfg topResidual subchunks
 where longrange = fold $ Arr.zipWith toUniformSampledLike
              (subdivideHomogenSampled nSubDivs lowpassed) modChunks
       lowpassed = homogenSampled (0,1) (simpleIIRLowpass lrBandwidth allYs)
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
 where HomogenSampled _ loRes _ = resampleHomogen lowpassed (0,1) nChunkMax
       transformedLr = invFourierTrafo loRes
       maxPosAmplitudeLr = Arr.foldl' maxℂ 0 transformedLr
       maxNegAmplitudeLr = Arr.foldl' minℂ 0 transformedLr
       μLr = maximum [maxPosAmplitudeLr, -maxNegAmplitudeLr, noiseLvl]
               / maxAllowedVal
       powerSpectrumLr = simpleIIRHighpass (fromIntegral nTot / subdivFreq)
                                $ Arr.map abs²ℂ transformedLr
       quantisedLr :: UArr.Vector (Complex c)
       quantisedLr = Arr.map (roundℂ . (/realToFrac μLr))
            $ filterFirstNPositivesOf powerSpectrumLr transformedLr
       subResiduals = transposeV
               $ map (subdivideHomogenSampled nSubDivs) residuals
       subchunks
        | nTot < localInfo
                     = Arr.empty
        | otherwise  = Arr.map (chunkFromUniform cfg) subResiduals
       maxAllowedVal = fromIntegral (maxBound :: c) / 8
       nTot = dynDimension lowpassed
       nSingleChunk = fromIntegral nTot / subdivFreq
       nFullChunk = ceiling $ fromIntegral nTot / fromIntegral subdivisionsSizeFactor
       nTaper = round $ fromIntegral nTot * subdivsOverlap
                             / fromIntegral subdivisionsSizeFactor
       taperStart, taperEnd :: UArr.Vector Double -> UArr.Vector Double
       taperStart = Arr.zipWith (*) $ Arr.generate nFullChunk
                      (\i -> if i<nTaper
                              then let x = fromIntegral i / fromIntegral nTaper
                                   in x^2 * (3 - 2*x)
                              else 1 )
       taperEnd = Arr.zipWith (*) $ Arr.generate nFullChunk
                      (\i -> let i' = nFullChunk - i
                             in if i' < nTaper
                              then let x = fromIntegral i' / fromIntegral nTaper
                                   in x^2 * (3 - 2*x)
                              else 1 )
       filterFirstNPositivesOf :: (Num n, UArr.Storable n) => UArr.Vector Double
                                   -> UArr.Vector n -> UArr.Vector n
       filterFirstNPositivesOf spect d = (`Arr.unfoldr`(0,infoPerStage))
           $ \(i,capacity) -> if capacity>0 && i<nTot
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
