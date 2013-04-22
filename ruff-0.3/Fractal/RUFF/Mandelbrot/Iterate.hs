{-# LANGUAGE BangPatterns, DeriveDataTypeable #-}
{- |
Module      :  Fractal.RUFF.Mandelbrot.Iterate
Copyright   :  (c) Claude Heiland-Allen 2011
License     :  BSD3

Maintainer  :  claudiusmaximus@goto10.org
Stability   :  unstable
Portability :  portable

Generic (slow) functions to iterate points.
-}

module Fractal.RUFF.Mandelbrot.Iterate
  ( Mode (..)
  , Iterate (..)
  , Output (..)
  , initial
  , iterate
  , iterates
  ) where

import Data.Data (Data)
import Data.Typeable (Typeable)
import Prelude hiding (iterate)
import Fractal.RUFF.Types.Complex (Complex(..), phase)

-- | Iteration mode.
data Mode = Simple | EscapeTime | DistanceEstimate
  deriving (Read, Show, Eq, Ord, Enum, Bounded, Data, Typeable)

-- | Iteration state.
data Iterate u r
  = IterSimple{ itc, itz :: !(Complex r), iterUser :: !u }
  | IterEscapeTime{ itc, itz :: !(Complex r), itn :: !Int, iterUser :: !u }
  | IterDistanceEstimate{ itc, itz, itdz :: !(Complex r), itn :: !Int, iterUser :: !u }
  deriving (Read, Show, Eq, Data, Typeable)

-- | Iteration initial state.
initial :: Num r => Mode -> u -> Complex r -> Iterate u r
{-# INLINABLE initial #-}
initial Simple           u c = IterSimple
  { itc = c, itz = 0 :+ 0,                         iterUser = u }
initial EscapeTime       u c = IterEscapeTime
  { itc = c, itz = 0 :+ 0,                itn = 0, iterUser = u }
initial DistanceEstimate u c = IterDistanceEstimate
  { itc = c, itz = 0 :+ 0, itdz = 0 :+ 0, itn = 0, iterUser = u }

-- | Iteration output.
data Output u r
  = OutSimple{ outUser :: !u }
  | OutEscapeTime{ escapeTime, finalAngle :: !r, outUser :: !u }
  | OutDistanceEstimate{ escapeTime, finalAngle, distanceEstimate :: !r, outUser :: !u }
  deriving (Read, Show, Eq, Ord, Data, Typeable)

-- | Iteration engine.
iterate :: (Ord r, Floating r) => Int -> Iterate u r -> Either (Iterate u r) (Output u r)
{-# INLINABLE iterate #-}
iterate n i@(IterSimple{ itc = cr :+ ci, itz = z0, iterUser = u }) = go 0 z0
  where
    go !m !z@(zr :+ zi)
      | m < n = let !zrr = zr * zr
                    !zii = zi * zi
                    !zri = zr * zi
                    !e = zrr + zii > 4
                in  if e then Right (OutSimple{ outUser = u})
                         else go (m + 1) ((zrr - zii + cr) :+ (2 * zri + ci))
      | otherwise = Left (i{ itz = z })
iterate n i@(IterEscapeTime{ itc = cr :+ ci, itz = z0, itn = n0, iterUser = u }) = go 0 z0
  where
    er = 65536
    er2 = er * er
    log2 = log 2
    go !m !z@(zr :+ zi)
      | m < n = let !zrr = zr * zr
                    !zii = zi * zi
                    !zri = zr * zi
                    !zz = zrr + zii
                    !e = zz > er2
                in  if e then Right (OutEscapeTime
                                { escapeTime = fromIntegral (n0 + m) + (log (log er) - log (log zz / 2)) / log2
                                , finalAngle = phase z
                                , outUser = u})
                         else go (m + 1) ((zrr - zii + cr) :+ (2 * zri + ci))
      | otherwise = Left (i{ itz = z, itn = n0 + n })
iterate !n !i@(IterDistanceEstimate{ itc = cr :+ ci, itz = z0, itdz = dz0, itn = n0, iterUser = u }) = go 0 z0 dz0
  where
    er = 65536
    er2 = er * er
    log2 = log 2
    go !m !z@(zr :+ zi) !dz@(dzr :+ dzi)
      | m < n = let !zrr = zr * zr
                    !zii = zi * zi
                    !zri = zr * zi
                    !zz = zrr + zii
                    !e = zz > er2
                    !zdzr = zr * dzr - zi * dzi
                    !zdzi = zr * dzi + zi * dzr
                    !dzdz = dzr * dzr + dzi * dzi
                in  if e then Right (OutDistanceEstimate
                                { escapeTime = fromIntegral (n0 + m) + (log (log er) - log (log zz / 2)) / log2
                                , finalAngle = phase z
                                , distanceEstimate = log zz * sqrt (zz / dzdz)
                                , outUser = u})
                         else go (m + 1) ((zrr - zii + cr) :+ (2 * zri + ci)) ((2 * zdzr + 1) :+ (2 * zdzi))
      | otherwise = Left (i{ itz = z, itdz = dz, itn = n0 + n })

-- | Iterate over a list.
iterates :: (Functor m, Monad m, Ord r, Floating r) => Int -> [Iterate u r] -> (Output u r -> m ()) -> m [Iterate u r]
{-# INLINABLE iterates #-}
iterates _ []     _   = return []
iterates n (x:xs) out = case iterate n x of
  Right o -> out o   >>  iterates n xs out
  Left  y -> (y:) `fmap` iterates n xs out
