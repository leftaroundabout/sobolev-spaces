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
                       , evalUnitL2
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

import Data.Monoid ((<>))

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

data SampleMode = DiscreteSineTransform

type UABSample c = (Integral c, Bounded c, UArr.Storable c)

-- | A square-integrable, bandlimited function with compact support
--   on the unit interval. 
data UnitL2 c s = UnitL2 {
    unitL2AmplitudeFactor :: !s
  , unitL2ExtremeSampleVals :: !(c,c)
  , unitL2LoFreqSamples :: UArr.Vector c
  , unitL2LoFreqSampleMode :: !SampleMode
  , unitL2Subdivisions :: BArr.Vector (UnitL2 c s)
  }

-- | Efficient only when partially applied to data array
--   and then evaluated at multiple x-points 
lfCubicSplineEval :: UArr.Vector Double -> Double -> Double
lfCubicSplineEval ys
  | n > 1   = lcse
  | n < 1   = const 0
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

sineTrafo :: UArr.Vector Double -> UArr.Vector Double
sineTrafo = prepareTrafos FFT.dst1 (\n αs -> αs <> UArr.replicate (n - Arr.length αs) 0)

invSineTrafo :: UArr.Vector Double -> UArr.Vector Double
invSineTrafo = prepareTrafos FFT.idst1 cubicResample

prepareTrafos :: (UArr.Storable a, UArr.Storable b)
                    => FFT.Transform a b
                     -> (Int -> UArr.Vector a -> UArr.Vector a)
                     -> UArr.Vector a -> UArr.Vector b
prepareTrafos transform resize = \αs
   -> let nmin = Arr.length αs
      in case Map.lookupGT nmin plans of
          Just (n,plan) -> FFT.execute plan
                            $ if n > nmin
                               then resize n αs
                               else αs
 where plans = Map.fromList [ (n, FFT.plan transform n)
                            | p <- [3..10]
                            , υ <- [0,1]
                            , let n = 2^p + υ*2^(p-1) -- [8,12,16,24,32,48..]
                            ]

evalUnitL2 :: (Integral c, UArr.Storable c) => UnitL2 c Double -> Double -> Double
evalUnitL2 (UnitL2 μ _ lf DiscreteSineTransform subdivs) = evalAt
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
                  (sineTrafo $ UArr.map ((*μ) . fromIntegral) lf)
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

fromUniformSampled :: ∀ c . UABSample c
        => SigSampleConfig
        -> UArr.Vector Double -> UnitL2 c Double
fromUniformSampled cfg@(SigSampleConfig nChunkMax
                           infoPerStage localInfo lrBandwidth noiseLvl) ys
  | nTot < localInfo
  , transformed <- invSineTrafo ys
  , maxPosAmplitude <- Arr.maximum transformed
  , maxNegAmplitude <- Arr.minimum transformed
  , μ <- maximum [maxPosAmplitude, -maxNegAmplitude, noiseLvl] / maxAllowedVal
     = UnitL2 μ (round $ maxNegAmplitude/μ, round $ maxPosAmplitude/μ)
              (Arr.map (round . (recip μ*)) transformed) DiscreteSineTransform
              Arr.empty
  | loRes <- cubicResample nChunkMax $ simpleIIRLowpass lrBandwidth ys
  , transformed <- invSineTrafo loRes
  , maxPosAmplitude <- Arr.maximum transformed
  , maxNegAmplitude <- Arr.minimum transformed
  , μ <- maximum [maxPosAmplitude, -maxNegAmplitude, noiseLvl] / maxAllowedVal
  , longrangePowerSpectrum <- simpleIIRHighpass (fromIntegral nTot / subdivFreq)
                                $ Arr.map (^2) transformed
  , loResQuantised <- Arr.map (round . (/μ))
            $ filterFirstNPositivesOf longrangePowerSpectrum transformed
  , backTransformed <- cubicResample nTot . sineTrafo
                         $ Arr.map ((*μ) . fromIntegral) loResQuantised
  , residual <- Arr.zipWith (-) ys backTransformed
  , chunks <- Arr.generate nSubDivs
               (\j -> let i = round (fromIntegral j * nSingleChunk)
                          i' = min (nTot-1)
                                $ round (fromIntegral (j+1) * nSingleChunk) + nTaper
                      in (if False && j>0 then taperStart else id)
                         . (if False && j<nSubDivs-1 then taperEnd else id)
                         $ Arr.slice i (i'-i) residual )
     = UnitL2 μ (round $ maxNegAmplitude/μ, round $ maxPosAmplitude/μ)
              loResQuantised DiscreteSineTransform
              (fromUniformSampled cfg <$> chunks)
 where maxAllowedVal = fromIntegral (maxBound :: c) / 8
       nTot = Arr.length ys
       nSingleChunk = fromIntegral nTot / subdivFreq
       nTaper = round $ fromIntegral nTot * subdivsOverlap
                             / fromIntegral subdivisionsSizeFactor
       taperStart, taperEnd :: UArr.Vector Double -> UArr.Vector Double
       taperStart = Arr.zipWith (*) $ Arr.generate (ceiling nSingleChunk)
                      (\i -> if i<nTaper
                              then let x = fromIntegral i / fromIntegral nTaper
                                   in x^2 * (3 - 2*x)
                              else 1 )
       taperEnd = Arr.zipWith (*) $ Arr.generate (ceiling nSingleChunk)
                      (\i -> let i' = ceiling nSingleChunk - i
                             in if i' < nTaper
                              then let x = fromIntegral i / fromIntegral nTaper
                                   in x^2 * (3 - 2*x)
                              else 1 )
       filterFirstNPositivesOf ::
            UArr.Vector Double -> UArr.Vector Double -> UArr.Vector Double
       filterFirstNPositivesOf spect d = (`Arr.unfoldr`(0,infoPerStage))
           $ \(i,capacity) -> if capacity>0 && i<nTot
                               then Just $ if spect!i > 0
                                            then (d!i, (i+1, capacity-1))
                                            else (0, (i+1, capacity))
                               else Nothing


