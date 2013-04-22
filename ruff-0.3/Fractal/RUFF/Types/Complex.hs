{-# LANGUAGE DeriveDataTypeable #-}
{- |
Module      :  Fractal.RUFF.Types.Complex
Copyright   :  (c) Claude Heiland-Allen 2011
License     :  BSD3

Maintainer  :  claudiusmaximus@goto10.org
Stability   :  unstable
Portability :  portable

Complex numbers without the 'RealFloat' constraint.
-}

module Fractal.RUFF.Types.Complex
  ( Complex((:+)), cis, mkPolar
  , realPart, imagPart, conjugate
  , magnitude2, magnitude, phase, polar
  ) where

import Data.Data (Data)
import Data.Typeable (Typeable)
import Data.Vec (NearZero(nearZero))

-- | Complex number type without the 'RealFloat' constraint.
data Complex r = !r :+ !r
  deriving (Read, Show, Eq, Data, Typeable)

instance Num r => Num (Complex r) where
  (x :+ y) + (u :+ v) = (x + u) :+ (y + v)
  (x :+ y) - (u :+ v) = (x - u) :+ (y - v)
  (x :+ y) * (u :+ v) = (x * u - y * v) :+ (x * v + y * u)
  negate (x :+ y) = negate x :+ negate y
  abs = error "Fractal.Types.Complex.Num.abs"
  signum = error "Fractal.Types.Complex.Num.signum"
  fromInteger n = fromInteger n :+ 0

instance Fractional r => Fractional (Complex r) where
  (x :+ y) / (u :+ v) = ((x * u + y * v) / d) :+ ((y * u - x * v) / d) where d = u * u + v * v
  fromRational r = fromRational r :+ 0

instance (Ord r, Floating r) => Floating (Complex r) where
  pi = pi :+ 0
  exp (x :+ y) = mkPolar (exp x) y
  log z = let (r, t) = polar z in log r :+ t
  sin (x :+ y) = (sin x * cosh y) :+        (cos x * sinh y)
  cos (x :+ y) = (cos x * cosh y) :+ negate (sin x * sinh y)
  tan z = sin z / cos z
  asin z = negate i * log (i * z + sqrt (1 - z*z)) where i = 0:+1
  acos z = negate i * log (z + sqrt (z*z - 1)) where i = 0:+1
  atan z = 1/2 * i * log ((1 - iz)/(1 + iz)) where i = 0:+1 ; iz = i * z
  sinh z = (exp z - exp (-z)) / 2
  cosh z = (exp z + exp (-z)) / 2
  tanh z = let ez2 = exp (2 * z) in (ez2 - 1) / (ez2 + 1)
  asinh z = log (z + sqrt (z*z + 1))
  acosh z = log (z + sqrt (z*z - 1))
  atanh z = 1/2 * log ((1 + z) / (1 - z))

instance NearZero r => NearZero (Complex r) where
  nearZero (r :+ i) = nearZero r && nearZero i

-- | Extract the real part.
realPart :: Complex r -> r
realPart (r :+ _) = r

-- | Extract the imaginary part.
imagPart :: Complex r -> r
imagPart (_ :+ i) = i

-- | Complex conjugate.
conjugate :: Num r => Complex r -> Complex r
conjugate (r :+ i) = r :+ negate i

-- | Complex magnitude squared.
magnitude2 :: Num r => Complex r -> r
magnitude2 (r :+ i) = r * r + i * i

-- | Complex magnitude.
magnitude :: Floating r => Complex r -> r
magnitude = sqrt . magnitude2

-- | Complex phase.
phase :: (Ord r, Floating r) => Complex r -> r
phase (r :+ i)
  | r > 0 && i > 0 =      atan (    i /     r)
  | r > 0 && i < 0 =    - atan (abs i /     r)
  | r < 0 && i > 0 = pi - atan (    i / abs r)
  | r < 0 && i < 0 =      atan (abs i / abs r) - pi
  | i > 0          =      pi / 2
  | i < 0          =    - pi / 2
  | r < 0          =      pi
  | otherwise      =      0

-- | Complex number with the given magnitude and phase.
mkPolar :: Floating r => r -> r -> Complex r
mkPolar r t = (r * cos t) :+ (r * sin t)

-- | Complex number with magnitude 1 and the given phase.
cis :: Floating r => r -> Complex r
cis t = cos t :+ sin t

-- | Convert to polar form.
polar :: (Ord r, Floating r) => Complex r -> (r, r)
polar z = (magnitude z, phase z)
