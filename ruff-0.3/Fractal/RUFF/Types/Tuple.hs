{-# LANGUAGE DeriveDataTypeable #-}
{- |
Module      :  Fractal.RUFF.Types.Tuple
Copyright   :  (c) Claude Heiland-Allen 2011
License     :  BSD3

Maintainer  :  claudiusmaximus@goto10.org
Stability   :  unstable
Portability :  portable

Strict tuples.
-}

module Fractal.RUFF.Types.Tuple
  ( Tuple2(..)
  ) where

import Data.Data (Data)
import Data.Typeable (Typeable)

-- | Strict 'Tuple2' type.
data Tuple2 l r = Tuple2 !l !r
  deriving (Read, Show, Eq, Ord, Data, Typeable)
