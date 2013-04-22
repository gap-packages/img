{-# LANGUAGE BangPatterns #-}
{- |
Module      :  Fractal.RUFF.Mandelbrot.Ray
Copyright   :  (c) Claude Heiland-Allen 2011
License     :  BSD3

Maintainer  :  claudiusmaximus@goto10.org
Stability   :  unstable
Portability :  portable

External angles define external rays which can be traced back from
the circle at infinity to parameters near the boundary of the Mandelbrot
Set.  Conversely, parameters near the boundary of the Mandelbrot Set can
be traced outwards to compute external angles.
-}
module Fractal.RUFF.Mandelbrot.Ray (externalRay, externalRayOut) where

import Data.Maybe (fromMaybe)

import Fractal.RUFF.Types.Complex (Complex, magnitude2, magnitude, phase, mkPolar)
import Fractal.RUFF.Mandelbrot.Address (Angle, double)

-- | Compute the external ray for an external angle with a given
--   accuracy, sharpness and starting radius.  For example:
--
-- > externalRay 1e-10 8 (2**24) (1/3)
--
--   The algorithm is based on Tomoki Kawahira's paper
--   /An algorithm to draw external rays of the Mandelbrot set/
--   <http://www.math.nagoya-u.ac.jp/~kawahira/programs/mandel-exray.pdf>.
--
externalRay :: (Ord r, Floating r) => r {- ^ accuracy -} -> Int {- ^ sharpness -} -> r {- ^ radius -} -> Angle {- ^ external angle -} -> [Complex r]
externalRay accuracy sharpness radius angle = map fst3 . iterate step $ (mkPolar radius (2 * pi * fromRational angle), accuracy * radius, (0, 0))
  where
    fst3 (x, _, _) = x
    -- step :: (NearZero r, Floating r) => (Complex r, (Int, Int)) -> (Complex r, (Int, Int))
    step (!c, !epsilon, (!k0, !j0))
      | j > sharpness = step (c, epsilon, (k0 + 1, 0))
      | otherwise =
          let c' = n c
              epsilon' = accuracy * magnitude (c' - c)
          in  (c', epsilon', (k0, j0 + 1))
      where
        epsilon2 = epsilon * epsilon
        k = k0 + 1
        j = j0 + 1
        m = (k - 1) * sharpness + j
        r = radius ** ((1/2) ** (fromIntegral m / fromIntegral sharpness))
        t = mkPolar (r ** (2 ** fromIntegral k0)) (2 * pi * fromRational (iterate double angle !! k0))
        n !z = let d = (cc - t) / dd in if not (magnitude2 d > epsilon2) then z else n (z - d)
          where
            (cc, dd) = ncnd k
            ncnd 1 = (z, 1)
            ncnd i = let (!nc, !nd) = ncnd (i - 1) in (nc * nc + z, 2 * nc * nd + 1)

-- | Compute the external ray outwards from a given parameter value.
--   If the result @rs@ satisfies:
--
--   > c = last rs
--   > magnitude c > radius
--
--   then the external angle is given by @t@:
--
--   > a = phase c / (2 * pi)
--   > t = a - fromIntegral (floor a)
--
externalRayOut :: (Ord r, Floating r, RealFrac r)
      => Int       {- ^ iterations -}
      -> r         {- ^ epsilon -}
      -> r         {- ^ accuracy -}
      -> Int       {- ^ sharpness -}
      -> r         {- ^ radius -}
      -> Complex r {- ^ parameter -}
      -> [Complex r]
externalRayOut maxIters epsilon accuracy sharpness radius = go (epsilon * epsilon)
  where
    radius2 = radius * radius
    iter !c !n !z
      | magnitude2 z > radius2 = Just (n, z)
      | n > maxIters = Nothing
      | otherwise = iter c (n + 1) (z * z + c)
    iterd !c !z !dz !m
      | m == 0 = (z, dz)
      | otherwise = iterd c (z * z + c) (2 * z * dz + 1) (m - 1)
    go !epsilon2 !c = (c :) . fromMaybe [] $ do
      (n, z) <- iter c 0 0
      let d = fromIntegral n - logBase 2 (log (magnitude2 z) / log radius2)
          d' = d - 1 / fromIntegral sharpness
          m = ceiling d'
          r = radius ** (2 ** (fromIntegral m - d'))
          a = phase z / (2 * pi)
          t = a - fromIntegral (floor a :: Int)
          k0 = mkPolar r (phase z)
          k1 = mkPolar r (pi *  t     )
          k2 = mkPolar r (pi * (t + 1))
          step !k !c0 = let (f, df) = iterd c0 0 0 m
                            dc = (f - k) / df
                            c0' = c0 - dc
                        in  c0 : if not (magnitude2 dc > epsilon2) then [] else step k c0'
          steps k = step k c
      if m == n
        then do
          return $ let c' = last $ steps k0 in go (accuracy * magnitude2 (c' - c)) c'
        else if m > 0 then do
          let (c1, c2) = last (steps k1 `zip` steps k2)
          (n1, _) <- iter c1 0 0
          (n2, _) <- iter c2 0 0
          let (c', n')
                | magnitude2 (c1 - c) < magnitude2 (c2 - c) = (c1, n1)
                | otherwise                                 = (c2, n2)
          return $ if n' == m then go (accuracy * magnitude2 (c' - c)) c' else []
        else return []
