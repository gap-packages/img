{-# LANGUAGE BangPatterns #-}
{- |
Module      :  Fractal.RUFF.Mandelbrot.Nucleus
Copyright   :  (c) Claude Heiland-Allen 2011
License     :  BSD3

Maintainer  :  claudiusmaximus@goto10.org
Stability   :  unstable
Portability :  portable

Mu-atom period, nucleus and bond point finding.
-}
module Fractal.RUFF.Mandelbrot.Nucleus (findPeriod, findNucleus, findBond, findInternal) where

import Data.List (genericIndex)
import Data.Maybe (listToMaybe)
import Fractal.RUFF.Types.Complex (Complex((:+)), mkPolar, magnitude2)

-- | Given the period and approximate location, successively refine
--   this estimate to a nucleus.
--
--   The algorithm is based on Robert Munafo's page
--   /Newton-Raphson method/
--   <http://mrob.com/pub/muency/newtonraphsonmethod.html>.
--
findNucleus :: (Floating r, Fractional r) => Integer {- ^ period -} -> Complex r {- ^ estimate -} -> [Complex r]
findNucleus p g = iterate go g
  where
    go !c =
      let step (!z, !d) = (z * z + c, 2 * z * d + 1)
          (zn, dn) = iterate step (0, 0) `genericIndex` p
      in  c - zn / dn

-- | Given the period and nucleus, find succesive refinements to the
--   bond point at a given internal angle.
--
--   The algorithm is based on ideas from
--   <http://mrob.com/pub/muency/derivative.html>.
--
findBond :: (Floating r, Fractional r) => Integer {- ^ period -} -> Complex r {- ^ nucleus -} -> r {- ^ angle -} -> [Complex r]
findBond p c0 a0 = findInternal p c0 1 a0

-- | Given the period and nucleus, find an interior point at a given internal
--   angle and radius in (0,1].
--
findInternal :: (Floating r, Fractional r) => Integer {- ^ period -} -> Complex r {- ^ nucleus -} -> r {- ^ radius -} -> r {- ^ angle -} -> [Complex r]
findInternal p c0 r0 a0 = snd `map` iterate go (c0, c0)
  where
    b0 = mkPolar r0 (2 * pi * a0)
    go (!z1, !c1) =
      let step (!a, !b, !c, !d, !e) =
              ( a * a + c1
              , 2 * a * b
              , 2 * (b * b + a * c)
              , 2 * a * d + 1
              , 2 * (a * e + b * d)
              )
          (an, bn, cn, dn, en) = iterate step (z1, 1, 0, 0, 0) `genericIndex` p
          y0 = z1 - an
          y1 = b0 - bn
          bn1 = bn - 1
          m = bn1 * en - dn * cn
          d0 = (y0 * en - dn * y1) / m
          d1 = (bn1 * y1 - y0 * cn) / m
      in  (z1 + d0, c1 + d1)

-- | Find the period of the lowest period nucleus inside a square.
--
--   The algorithm is based on Robert Munafo's page,
--   /Finding the Period of a mu-Atom/
--   <http://mrob.com/pub/muency/period.html>.
--
findPeriod :: (Floating r, Ord r) => Integer {- ^ maximum period -} -> r {- ^ radius -} -> Complex r {- ^ center -} -> Maybe Integer
findPeriod m r c =
  let cs = [ c + (r:+r), c + (r:+(-r)), c + ((-r):+(-r)), c + ((-r):+r) ]
      zs = iterate (zipWith (\cc z -> z * z + cc) cs) [0,0,0,0]
      -- kludge = if r > 0 then 1 else 2 -- fixes space leak from long literal list (CAF?)
  in  fmap fst . listToMaybe . dropWhile (not . straddlesOrigin . snd) . takeWhile (all ((< 65536) . magnitude2) . snd) . zip [{-kludge + 0 - kludge-} 0 .. m ] $ zs

straddlesOrigin :: (Ord r, Num r) => [Complex r] -> Bool
straddlesOrigin ps = odd . length . filter id . zipWith positiveReal ps $ (drop 1 ps ++ take 1 ps)

positiveReal :: (Ord r, Num r) => Complex r -> Complex r -> Bool
positiveReal (u:+v) (x:+y)
  | v < 0 && y < 0 = False
  | v > 0 && y > 0 = False
  | (u * (y - v) - v * (x - u)) * (y - v) > 0 = True
  | otherwise = False
