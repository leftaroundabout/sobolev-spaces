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

module Math.FunctionalAnalysis.L2Function.R1 (
                         UnitL2
            ) where

import qualified Data.Vector.Unboxed as UArr
import qualified Data.Vector as BArr
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as Arr
-- import Numeric.FFT.Vector.Invertible


loResSampleCount :: Int
loResSampleCount = 64

subdivisionsSizeFactor :: Int
subdivisionsSizeFactor = 8

nSubDivs :: Int
nSubDivs = 10

data SampleMode = CubicInterpolation

-- | A square-integrable, bandlimited function with compact support
--   on the unit interval. 
data UnitL2 c s = UnitL2 {
    unitL2AmplitudeFactor :: !s
  , unitL2MinimumSampleVal, unitL2MaximumSampleVal :: !c
  , unitL2LoFreqSamples :: UArr.Vector c
  , unitL2LoFreqSampleMode :: !SampleMode
  , unitL2Subdivisions :: BArr.Vector (UnitL2 c s)
  }

-- | Efficient only when partially applied to data array
--   and then evaluated at multiple x-points 
lfCubicSplineEval :: (Integral c, UArr.Unbox c)
            => (Double, UArr.Vector c) -> Double -> Double
lfCubicSplineEval (μ,ys)
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
              y₀ = μ * fromIntegral (ys!(i-1))
              y₁ = μ * fromIntegral (ys!i)
              ð₀ = μ/2 * δys!(i-1)
              ð₁ = μ/2 * δys!i
       δys = UArr.imap (\i _
                  -> if i<1
                      then fromIntegral (ys!1)
                      else if i<n-1
                            then fromIntegral (ys!(i+1)) - fromIntegral (ys!(i-1))
                            else -fromIntegral (ys!(n-2))
              ) ys
       n = Arr.length ys
       δx = 1 / fromIntegral loResSampleCount

evalUnitL2 :: (Integral c, Bounded c) => UnitL2 c Double -> Double -> Double
evalUnitL2 = undefined

