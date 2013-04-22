{-# LANGUAGE BangPatterns, DeriveDataTypeable #-}
{- |
Module      :  Fractal.RUFF.Mandelbrot.Image
Copyright   :  (c) Claude Heiland-Allen 2011
License     :  BSD3

Maintainer  :  claudiusmaximus@goto10.org
Stability   :  unstable
Portability :  portable

Generic (slow) functions to render images.

-}

module Fractal.RUFF.Mandelbrot.Image
  ( simpleImage, complexImage, imageLoop, coordinates, ascii, unicode
  , Channel(..), Coordinates, border
  ) where

import Control.Monad.ST (ST)
import Data.Array.ST (newArray, writeArray, runSTUArray)
import Data.STRef (STRef, newSTRef, readSTRef, writeSTRef)
import Data.Array.Unboxed (UArray, (!), bounds, range, amap, ixmap)

import Data.Ix (Ix)
import Data.Data (Data)
import Data.Typeable (Typeable)

import Fractal.RUFF.Types.Complex (Complex((:+)), magnitude)
import Fractal.RUFF.Types.Tuple (Tuple2(Tuple2))
import Fractal.RUFF.Mandelbrot.Iterate (iterates, initial, Mode(Simple, DistanceEstimate), Iterate(), Output(OutSimple, OutDistanceEstimate), escapeTime, distanceEstimate, finalAngle, outUser)

-- | Render an image with the 'Simple' algorithm.  The iteration count is
--   doubled until the image is good enough, or the fixed maximum iteration
--   count is reached.
--
-- > putStr . unicode $ simpleImage (coordinates 100 100 ((-1.861):+0) (0.001)) 1000000000
simpleImage :: (Ord r, Floating r) => Coordinates r {- ^ coordinates -} -> Int {- ^ max iterations -} -> UArray (Int, Int) Bool {- ^ image -}
{-# INLINABLE simpleImage #-}
simpleImage (bs, cs) n0 = runSTUArray $ do
    a <- newArray bs True
    s <- newSTRef (0 :: Int)
    imageLoop s a n0 0 False 64 i0s (out s a)
  where
    i0s = map (uncurry $ initial Simple) cs
    out s a (OutSimple{ outUser = Tuple2 j i }) = do
      writeArray a (j, i) False
      modifySTRef' s (+ 1)
    out _ _ _ = return ()
 
-- | Render an image with the 'DistanceEstimate' algorithm.  The iteration count is
--   doubled until the image is good enough, or the fixed maximum iteration
--   count is reached.  The output values are converted to 'Float'.
--
-- > putStr . unicode . border $ complexImage (coordinates 100 100 ((-1.861):+0) (0.001)) 1000000000
complexImage :: (Ord r, Real r, Floating r) => Coordinates r {-^ coordinates -} -> Int {- ^ max iterations -} -> UArray (Int, Int, Channel) Float {- ^ image -}
{-# INLINABLE complexImage #-}
complexImage (((jlo,ilo),(jhi,ihi)), cs) !n0 = runSTUArray $ do
    a <- newArray bs (-1)
    s <- newSTRef (0 :: Int)
    imageLoop s a n0 0 False 64 i0s (out s a)
  where
    bs = ((jlo,ilo,minBound), (jhi,ihi,maxBound))
    (_, cx0):(_, cx1):_ = cs
    pixelSpacing = magnitude (cx1 - cx0)
    i0s = map (uncurry $ initial DistanceEstimate) cs
    out !s !a (OutDistanceEstimate{ escapeTime = et, distanceEstimate = de, finalAngle = fa, outUser = Tuple2 j i }) = {-# SCC "complexImage.out" #-} do
      writeArray a (j, i, EscapeTime) (realToFrac et)
      writeArray a (j, i, DistanceEstimate') (realToFrac (de / pixelSpacing))
      writeArray a (j, i, FinalAngle) (realToFrac fa)
      modifySTRef' s (+ 1)
    out _ _ _ = return ()

-- | Channels in an image.
data Channel = EscapeTime {- ^ continuous dwell -} | DistanceEstimate' {- ^ normalized to pixel spacing -} | FinalAngle {- ^ in [-pi,pi] -}
  deriving (Eq, Ord, Enum, Bounded, Ix, Read, Show, Data, Typeable)

-- | Image rendering loop.
imageLoop :: (Ord r, Floating r) => STRef s Int {- ^ escapees -} -> a {- ^ output array -} -> Int {- ^ max iterations -} -> Int {- ^ iterations -} -> Bool {- ^ prior escapees -} -> Int {- ^ iterations this phase -} -> [Iterate u r] {- ^ iterates -} -> (Output u r -> ST s ()) {- ^ output callback -} -> ST s a {- ^ output array as given -}
{-# INLINABLE imageLoop #-}
imageLoop s a !n0 !n1 !f1 !m1 is1 out = loop f1 n1 m1 is1
  where
    loop !f !n !m is = do
      writeSTRef s 0
      is' <- iterates m is out
      o <- readSTRef s
      if null is || (f && o == 0) || n > n0 then return a else loop (f || o > 0) (n + m) (m * 2) is'

-- | Image bounds and coordinates.
type Coordinates r = (((Int,Int),(Int,Int)), [(Tuple2 Int Int, Complex r)])

-- | The parameter plane coordinates for an image, with bounds.
coordinates :: (Ord r, Floating r) => Int {- ^ width -} -> Int {- ^ height -} -> Complex r {- ^ center -} -> r {- ^ size -} -> Coordinates r
{-# INLINABLE coordinates #-}
coordinates !width !height !(c0r :+ c0i) !r0 = (bs, cs)
  where
    bs = ((0, 0), (height - 1, width - 1))
    cs =  [ (Tuple2 j i, c)
          | (j,i) <- range bs
          , let y = (fromIntegral j - h) / h
          , let x = (fromIntegral i - w) / h
          , let ci = c0i + r0 * y
          , let cr = c0r + r0 * x
          , let c = cr :+ ci
          ]
    w = fromIntegral $ width  `div` 2
    h = fromIntegral $ height `div` 2

-- | Convert a distance estimate image to a near-boundary bit array.
--   The input image must have a DistanceEstimate' channel.
border :: UArray (Int, Int, Channel) Float {- ^ image -} -> UArray (Int, Int) Bool
border a = amap (\x -> x > 0 && x < 1) . ixmap bs (\(j, i) -> (j, i, DistanceEstimate')) $ a
  where
    ((jlo, ilo, _), (jhi, ihi, _)) = bounds a
    bs = ((jlo, ilo), (jhi, ihi))

-- | Convert a bit array to ascii graphics.
ascii :: UArray (Int, Int) Bool {- ^ image -} -> String {- ^ ascii -}
ascii a = unlines . map concat $ [ [ b (a ! (j, i)) | i <- [ ilo .. ihi ] ] | j <- [ jhi, jhi - 1 .. jlo ] ]
  where
    ((jlo, ilo), (jhi, ihi)) = bounds a
    b False = "  "
    b True  = "##"

-- | Convert a bit array to unicode block graphics.
unicode :: UArray (Int, Int) Bool {- ^ image -} -> String {- ^ unicode -}
unicode a = unlines [ [ b (a ! (j, i)) (a ! (j - 1, i)) | i <- [ ilo .. ihi ] ] | j <- [ jhi, jhi - 2 .. jlo ] ]
  where
    ((jlo, ilo), (jhi, ihi)) = bounds a
    b False False = ' '
    b True False = '\x2580'
    b False True = '\x2584'
    b True True = '\x2588'

-- | Strict version of 'modifySTRef'.
modifySTRef' :: STRef s a -> (a -> a) -> ST s ()
{-# INLINABLE modifySTRef' #-}
modifySTRef' s f = do
  x <- readSTRef s
  writeSTRef s $! f x
