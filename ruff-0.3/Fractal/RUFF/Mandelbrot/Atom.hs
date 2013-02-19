{-# LANGUAGE BangPatterns #-}
{- |
Module      :  Fractal.RUFF.Mandelbrot.Atom
Copyright   :  (c) Claude Heiland-Allen 2011
License     :  BSD3

Maintainer  :  claudiusmaximus@goto10.org
Stability   :  unstable
Portability :  portable

Mu-atom coordinate and address algorithms.
-}
module Fractal.RUFF.Mandelbrot.Atom
  ( MuAtom(..)
  , FindAtom(..), findAtom, findAtom', findAtom_
  , FindAddress(..), findAddress, findAddress', findAddress_
  , Locate(..), locate, locate', locate_
  ) where

import Control.Arrow ((***))
import Data.Maybe (listToMaybe)
import Data.Ratio ((%))
import Data.Vec (NearZero, nearZero)

import Fractal.RUFF.Mandelbrot.Address (AngledInternalAddress, Angle, splitAddress, addressPeriod, externalAngles, angledInternalAddress)
import Fractal.RUFF.Mandelbrot.Nucleus (findNucleus, findBond, findPeriod)
import Fractal.RUFF.Mandelbrot.Ray (externalRay, externalRayOut)
import Fractal.RUFF.Types.Complex (Complex, magnitude, magnitude2, phase, mkPolar)

-- | Mu-atom properties.
data MuAtom r = MuAtom
  { muNucleus :: !(Complex r)
  , muSize    :: !Double
  , muOrient  :: !Double
  , muPeriod  :: !Integer
  }
  deriving (Read, Show, Eq)

-- | Progress updates for 'findAtom'.
data FindAtom r
  = AtomSplitTodo
  | AtomSplitDone AngledInternalAddress [Angle]
  | AtomAnglesTodo
  | AtomAnglesDone !Rational !Rational
  | AtomRayTodo
  | AtomRay !Integer
  | AtomRayDone !(Complex r)
  | AtomNucleusTodo
  | AtomNucleus !Integer
  | AtomNucleusDone !(Complex r)
  | AtomBondTodo
  | AtomBond !Integer
  | AtomBondDone !(Complex r)
  | AtomSuccess !(MuAtom r)
  | AtomFailed
  deriving (Read, Show, Eq)

isAtomSuccess :: FindAtom r -> Bool
isAtomSuccess (AtomSuccess _) = True
isAtomSuccess _ = False

fromAtomSuccess :: FindAtom r -> Maybe (MuAtom r)
fromAtomSuccess (AtomSuccess m) = Just m
fromAtomSuccess _ = Nothing

-- | Try to find an atom, providing progress updates.
findAtom :: (Floating r, NearZero r, Real r) => AngledInternalAddress -> [FindAtom r]
findAtom addr = AtomSplitTodo :
  let (!iaddr, !caddr) = splitAddress addr
      !p = addressPeriod iaddr
  in  AtomSplitDone iaddr caddr : AtomAnglesTodo : case externalAngles iaddr of
    Nothing -> [AtomFailed]
    Just (!lo, !hi) -> AtomAnglesDone lo hi : AtomRayTodo :
      let sharpness = 8
          er = 65536
          accuracy = 1e-10
          ok w = magnitude2 w < 2 * er ^ (2::Int) -- NaN -> False
          rayl = externalRay accuracy sharpness er lo
          rayh = externalRay accuracy sharpness er hi
          ray' = takeWhile (uncurry (&&) . (ok *** ok) . snd) $ [ 1 .. ] `zip` (rayl `zip` rayh)
          rgo []  _ = [AtomFailed]
          rgo [_] _ = [AtomFailed]
          rgo ((i, (xl, xh)):m@((_, (yl, yh)):_)) f
            | i < fromIntegral sharpness * (p + 16) = AtomRay i : rgo m f
            | dl + dh > dx + dy = AtomRay i : rgo m f
            | otherwise = AtomRayDone x : f x
            where
              x  = 0.5 * (xl + xh)
              dl = magnitude2 (xl - yl)
              dh = magnitude2 (xh - yh)
              dx = magnitude2 (xl - xh)
              dy = magnitude2 (yl - yh)
      in  rgo ray' $ \rayend -> AtomNucleusTodo :
        let nuc = findNucleus p rayend
            nuc' = takeWhile (ok . snd) $ [ 1 .. ] `zip` nuc
            ngo []  _ = [AtomFailed]
            ngo [_] _ = [AtomFailed]
            ngo ((i, x):m@((_, y):_)) f
              | not (nearZero (x - y)) = AtomNucleus i : ngo m f
              | otherwise = AtomNucleusDone x : f x
        in  ngo nuc' $ \nucleus -> AtomBondTodo :
          let bnd = findBond p nucleus 0.5
              bnd' = takeWhile (ok . snd) $ [ 1.. ] `zip` bnd
              bgo []  _ = [AtomFailed]
              bgo [_] _ = [AtomFailed]
              bgo ((i, x):m@((_, y):_)) f
                | not (nearZero (x - y)) = AtomBond i : bgo m f
                | otherwise = AtomBondDone x : f x
          in  bgo bnd' $ \bond ->
            let delta  = bond - nucleus
                size   = realToFrac $ magnitude delta / 0.75 -- FIXME check for island-hood
                orient = realToFrac $ phase delta
                atom   = MuAtom{ muNucleus = nucleus, muSize = size, muOrient = orient, muPeriod = p }
            in  if 10 > size && size > 0 then [AtomSuccess atom] else [AtomFailed]

-- | Find the first success in the progress list.
findAtom' :: [FindAtom r] -> Maybe (MuAtom r)
findAtom' ps = fromAtomSuccess =<< listToMaybe (filter isAtomSuccess ps)

-- | Find an atom from its address.
findAtom_ :: (Floating r, NearZero r, Real r) => AngledInternalAddress -> Maybe (MuAtom r)
findAtom_ = findAtom' . findAtom

-- | Progress updates for 'findAddress'.
data FindAddress r
  = AddressCuspTodo
  | AddressCuspDone !(Complex r)
  | AddressDwellTodo
  | AddressDwell !Integer
  | AddressDwellDone !Integer
  | AddressRayOutTodo
  | AddressRayOut !Double
  | AddressRayOutDone !(Complex r)
  | AddressExternalTodo
  | AddressExternalDone !Rational
  | AddressAddressTodo
  | AddressSuccess AngledInternalAddress
  | AddressFailed
  deriving (Read, Show, Eq)

isAddressSuccess :: FindAddress r -> Bool
isAddressSuccess (AddressSuccess _) = True
isAddressSuccess _ = False

fromAddressSuccess :: FindAddress r -> Maybe AngledInternalAddress
fromAddressSuccess (AddressSuccess a) = Just a
fromAddressSuccess _ = Nothing

-- | Try to find an address, providing progress updates.
findAddress :: (Floating r, NearZero r, Real r, RealFrac r) => MuAtom r -> [FindAddress r]
findAddress mu = AddressCuspTodo :
  let cusp = muNucleus mu - mkPolar (realToFrac (muSize mu)) (realToFrac ((muOrient mu)))
      er = 65536
      er2 = er * er
  in  AddressCuspDone cusp : AddressDwellTodo :
    let dgo z n f = AddressDwell n : if magnitude2 z > er2 then f n else dgo (z * z + cusp) (n + 1) f
    in  dgo 0 0 $ \n -> AddressDwellDone n : AddressRayOutTodo :
      let rgo ((i,!_):izs@(_:_)) f = AddressRayOut (fromIntegral i / (fromIntegral sharpness * fromIntegral n)) : rgo izs f
          rgo [(_,!z)] f | magnitude2 z > er2 = AddressRayOutDone z : f z
          rgo _ _ = [AddressFailed]
          accuracy = 1e-16
          sharpness = 16
          epsilon0 = realToFrac (muSize mu) * accuracy
      in  rgo ([(1 :: Integer) ..] `zip` externalRayOut (fromIntegral n + 100) epsilon0 accuracy sharpness er cusp) $ \rend -> AddressExternalTodo :
        let den = 2 ^ muPeriod mu - 1
            num' = fromIntegral den * warp (phase rend / (2 * pi))
            num = round num'
            warp t
              | t > 0 = t
              | otherwise = t + 1
            angle = num % den
        in  AddressExternalDone angle : AddressAddressTodo : case angledInternalAddress angle of
              Nothing -> [AddressFailed]
              Just addr -> if addressPeriod addr /= muPeriod mu then [AddressFailed] else [AddressSuccess addr]

-- | Find the first success in the progress list.
findAddress' :: [FindAddress r] -> Maybe AngledInternalAddress
findAddress' ps = fromAddressSuccess =<< listToMaybe (filter isAddressSuccess ps)

-- | Find an address for a mu-atom.
findAddress_ :: (Floating r, NearZero r, Real r, RealFrac r) => MuAtom r -> Maybe AngledInternalAddress
findAddress_ = findAddress' . findAddress

-- | Progress updates for 'locate'.
data Locate r
  = LocateScanTodo
  | LocateScan
  | LocateScanDone !Integer
  | LocateNucleusTodo
  | LocateNucleus !Integer
  | LocateNucleusDone !(Complex r)
  | LocateBondTodo
  | LocateBond !Integer
  | LocateBondDone !(Complex r)
  | LocateSuccess !(MuAtom r)
  | LocateFailed
  deriving (Read, Show, Eq)

isLocateSuccess :: Locate r -> Bool
isLocateSuccess (LocateSuccess _) = True
isLocateSuccess _ = False

fromLocateSuccess :: Locate r -> Maybe (MuAtom r)
fromLocateSuccess (LocateSuccess m) = Just m
fromLocateSuccess _ = Nothing

-- | Try to find an atom close to a coordinate, providing progress updates.
locate :: (Floating r, NearZero r, Real r) => Complex r {- ^ center -} -> Double {- ^ radius -} -> [Locate r]
locate c r = LocateScanTodo : LocateScan : case findPeriod 10000000 (realToFrac r) c of  -- FIXME hardcoded
  Nothing -> [LocateFailed]
  Just p -> LocateScanDone p : LocateNucleusTodo :
    let ok w = magnitude2 w < 16 -- NaN -> False
        nuc = findNucleus p c
        nuc' = takeWhile (ok . snd) $ [ 1 .. ] `zip` nuc
        ngo []  _ = [LocateFailed]
        ngo [_] _ = [LocateFailed]
        ngo ((i, x):m@((_, y):_)) f
          | not (nearZero (x - y)) = LocateNucleus i : ngo m f
          | otherwise = LocateNucleusDone x : f x
    in  ngo nuc' $ \nucleus ->LocateBondTodo :
      let bnd = findBond p nucleus 0.5
          bnd' = takeWhile (ok . snd) $ [ 1.. ] `zip` bnd
          bgo []  _ = [LocateFailed]
          bgo [_] _ = [LocateFailed]
          bgo ((i, x):m@((_, y):_)) f
            | not (nearZero (x - y)) = LocateBond i : bgo m f
            | otherwise = LocateBondDone x : f x
      in  bgo bnd' $ \bond ->
        let delta  = bond - nucleus
            size   = realToFrac $ magnitude delta / 0.75 -- FIXME only valid for islands
            orient = realToFrac $ phase delta
            atom   = MuAtom{ muNucleus = nucleus, muSize = size, muOrient = orient, muPeriod = p }
        in  if 10 > size && size > 0 then [LocateSuccess atom] else [LocateFailed]

-- | Find the first success in the progress list.
locate' :: [Locate r] -> Maybe (MuAtom r)
locate' ps = fromLocateSuccess =<< listToMaybe (filter isLocateSuccess ps)

-- | Find an atom close to a coordinate.
locate_ :: (Floating r, NearZero r, Real r) => Complex r -> Double -> Maybe (MuAtom r)
locate_ c r = locate' (locate c r)
