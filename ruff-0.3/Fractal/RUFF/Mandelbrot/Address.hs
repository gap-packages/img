{-# LANGUAGE DeriveDataTypeable #-}
{- |
Module      :  Fractal.RUFF.Mandelbrot.Address
Copyright   :  (c) Claude Heiland-Allen 2010-2011
License     :  BSD3

Maintainer  :  claudiusmaximus@goto10.org
Stability   :  unstable
Portability :  portable

External angles give rise to kneading sequences under the angle doubling
map.  Internal addresses encode kneading sequences in human-readable form,
when extended to angled internal addresses they distinguish hyperbolic
components in a concise and meaningful way.

The algorithms are mostly based on Dierk Schleicher's paper
/Internal Addresses Of The Mandelbrot Set And Galois Groups Of Polynomials (version of February 5, 2008)/
<http://arxiv.org/abs/math/9411238v2>.
-}

module Fractal.RUFF.Mandelbrot.Address
  ( Angle, double, wrap, prettyAngle, prettyAngles
  , Knead(..), kneadChar
  , Kneading(..), prettyKneading, kneading, period, unwrap, associated, upper, lower
  , InternalAddress(..), prettyInternalAddress, internalAddress, internalFromList, internalToList
  , AngledInternalAddress(..), prettyAngledInternalAddress, angledInternalAddress, angledFromList, angledToList
  , externalAngles, stripAngles, splitAddress, joinAddress, addressPeriod
  , parseAngle, parseAngles, parseKnead, parseKneading, parseInternalAddress, parseAngledInternalAddress
  ) where

import Data.Data (Data())
import Data.Typeable (Typeable())
import Control.Monad (guard)
import Control.Monad.Identity (Identity())
import Data.Char (digitToInt)
import Data.Bits (testBit)
import Data.List (genericDrop, genericIndex, genericLength, genericReplicate, genericSplitAt, genericTake, foldl')
import Data.Maybe (isJust, listToMaybe)
import Data.Ratio ((%), numerator, denominator)
import Text.Parsec (ParsecT(), choice, digit, eof, many, many1, runP, sepEndBy, string, try)

-- | Angle as a fraction of a turn, usually in [0, 1).
type Angle = Rational

-- | Convert to human readable form.
prettyAngle :: Angle -> String
prettyAngle a = show (numerator a) ++ " / " ++ show (denominator a)

-- | Convert to human readable form.
prettyAngles :: [Angle] -> String
prettyAngles [] = ""
prettyAngles [a] = show (numerator a) ++ "/" ++ show (denominator a)
prettyAngles (a:as) = show (numerator a) ++ "/" ++ show (denominator a) ++ " " ++ prettyAngles as

-- | Wrap an angle into [0, 1).
wrap :: Angle -> Angle
wrap a
  | f < 0 = 1 + f
  | otherwise = f
  where
    (_, f) = properFraction a :: (Integer, Angle)

-- | Angle doubling map.
double :: Angle -> Angle
double a = wrap (2 * a)

-- | Binary representation of a (pre-)periodic angle.
type BinAngle = ([Bool], [Bool])

-- | Convert an angle from binary representation.
unbinary :: BinAngle -> Angle
unbinary (pre, per)
  | n == 0 = bits pre % (2 ^ m)
  | otherwise = (bits pre % (2 ^ m)) + (bits per % (2 ^ m * (2 ^ n - 1)))
  where
    m = length pre
    n = length per

-- | Convert a list of bits to an integer.
bits :: [Bool] -> Integer
bits = foldl' (\ a b -> 2 * a + if b then 1 else 0) 0

-- | Convert an angle to binary representation.
binary :: Angle -> BinAngle
binary a
  | a == 0 = ([], [])
  | even (denominator a) =
      let (pre, per) = binary (double a)
          b = a >= 1/2
      in  (b:pre, per)
  | otherwise =
      let (t, p) = head . dropWhile ((1 /=) . denominator . fst) . map (\q -> (a * (2^q - 1), q)) $ [ (1 :: Int) ..]
          s = numerator t
          n = fromIntegral p
          per = [ s `testBit` i | i <- [n - 1, n - 2 .. 0] ]
      in  ([], per)

-- | Tuning transformation for binary represented periodic angles.
--   Probably only valid for angle pairs presenting ray pairs.
btune :: BinAngle -> (BinAngle, BinAngle) -> BinAngle
btune (tpre, tper) (([], per0), ([], per1)) = (concatMap f tpre, concatMap f tper)
  where
    f False = per0
    f True  = per1
btune _ _ = error "btune: can't handle pre-periods"

-- | Tuning transformation for angles.
tune :: Angle -> (Angle, Angle) -> Angle
tune t (t0, t1) = unbinary $ btune (binary t) (binary t0, binary t1)

-- | Elements of kneading sequences.
data Knead
  = Zero
  | One
  | Star
  deriving (Read, Show, Eq, Ord, Enum, Bounded, Data, Typeable)

-- | Knead character representation.
kneadChar :: Knead -> Char
kneadChar Zero = '0'
kneadChar One  = '1'
kneadChar Star = '*'

-- | Kneading sequences.  Note that the 'Aperiodic' case has an infinite list.
data Kneading
  = Aperiodic [Knead]
  | PrePeriodic [Knead] [Knead]
  | StarPeriodic [Knead]
  | Periodic  [Knead]
  deriving (Read, Show, Eq, Ord, Data, Typeable)

-- | Kneading sequence as a string.  The 'Aperiodic' case is truncated arbitrarily.
prettyKneading :: Kneading -> String
prettyKneading (Aperiodic ks) = map kneadChar (take 17 ks) ++ "..."
prettyKneading (PrePeriodic us vs) = map kneadChar us ++ "(" ++ map kneadChar vs ++ ")"
prettyKneading (StarPeriodic vs) = "(" ++ map kneadChar vs ++ ")"
prettyKneading (Periodic vs) = "(" ++ map kneadChar vs ++ ")"

-- | The kneading sequence for an external angle.
kneading :: Angle -> Kneading
kneading a0'
  | a0 == 0 = StarPeriodic [Star]
  | otherwise = fst kneads
  where
    a0 = wrap a0'
    lo =  a0      / 2
    hi = (a0 + 1) / 2
    kneads = kneading' 1 (double a0)
    ks = (a0, One) : snd kneads
    kneading' :: Integer -> Angle -> (Kneading, [(Angle, Knead)])
    kneading' n a
      | isJust i = case i of
          Just 0 -> case last qs of
            Star -> (StarPeriodic qs, [])
            _    -> (Periodic qs, [])
          Just j -> let (p, q) = genericSplitAt j qs
                    in (PrePeriodic p q, [])
          _ -> error "Fractal.Mandelbrot.Address.kneading (isJust -> Nothing?)"
      | a == lo          = ((a, Star):) `mapP` k
      | a == hi          = ((a, Star):) `mapP` k
      | lo < a && a < hi = ((a, One ):) `mapP` k
      | hi < a || a < lo = ((a, Zero):) `mapP` k
      | otherwise = error "Fractal.Mandelbrot.Address.kneading (unmatched?)"
      where
        k = kneading' (n+1) (double a)
        ps = genericTake n ks
        qs = map snd ps
        i = fmap fst . listToMaybe . filter ((a ==) . fst . snd) . zip [(0 :: Integer) ..] $ ps
        mapP f ~(x, y) = (x, f y)

-- | The period of a kneading sequence, or 'Nothing' when it isn't periodic.
period :: Kneading -> Maybe Integer
period (StarPeriodic k) = Just (genericLength k)
period (Periodic k) = Just (genericLength k)
period _ = Nothing

rho :: Kneading -> Integer -> Integer
rho v r | r >= 1 && fmap (r`mod`) (period v) /= Just 0 = ((1 + r) +) . genericLength . takeWhile id . zipWith (==) vs . genericDrop r $ vs
        | otherwise = rho v (r + 1)
  where
    vs = unwrap v

-- | Unwrap a kneading sequence to an infinite list.
unwrap :: Kneading -> [Knead]
unwrap (Aperiodic vs) = vs
unwrap (PrePeriodic us vs) = us ++ cycle vs
unwrap (StarPeriodic vs) = cycle vs
unwrap (Periodic vs) = cycle vs

orbit :: (a -> a) -> a -> [a]
orbit = iterate

-- | Internal addresses are a non-empty sequence of strictly increasing
--   integers beginning with '1'.
data InternalAddress = InternalAddress [Integer]
  deriving (Read, Show, Eq, Ord, Data, Typeable)

-- | Internal address as a string.
prettyInternalAddress :: InternalAddress -> String
prettyInternalAddress (InternalAddress [])  = error "Fractal.Mandelbrot.Address.InternalAddress.pretty"
prettyInternalAddress (InternalAddress [x]) = show x
prettyInternalAddress (InternalAddress (x:ys)) = show x ++ " " ++ prettyInternalAddress (InternalAddress ys)

-- | Construct a valid 'InternalAddress', checking the precondition.
internalFromList :: [Integer] -> Maybe InternalAddress
internalFromList x0s@(1:_) = InternalAddress `fmap` fromList' 0 x0s
  where
    fromList' n [x]    | x > n = Just [x]
    fromList' n (x:xs) | x > n = (x:) `fmap` fromList' x xs
    fromList' _ _ = Nothing
internalFromList _ = Nothing

-- | Extract the sequence of integers.
internalToList :: InternalAddress -> [Integer]
internalToList (InternalAddress xs) = xs

-- | Construct an 'InternalAddress' from a kneading sequence.
internalAddress :: Kneading -> Maybe InternalAddress
internalAddress (StarPeriodic [Star])      = Just (InternalAddress [1])
internalAddress (StarPeriodic v@(One:_))   = Just . InternalAddress . address'per (genericLength v) $ v
internalAddress (Periodic     v@(One:_))   = Just . InternalAddress . address'per (genericLength v) $ v
internalAddress k@(Aperiodic    (One:_))   = Just . InternalAddress . address'inf . unwrap $ k
internalAddress _ = Nothing

address'inf :: [Knead] -> [Integer]
address'inf v = address' v

address'per :: Integer -> [Knead] -> [Integer]
address'per p v = takeWhile (<= p) $ address' v

address' :: [Knead] -> [Integer]
address' v = address'' 1 [One]
  where
    address'' sk vk = sk : address'' sk' vk'
      where
        sk' = (1 +) . genericLength . takeWhile id . zipWith (==) v . cycle $ vk
        vk' = genericTake sk' (cycle v)

-- | A star-periodic kneading sequence's upper and lower associated
--   kneading sequences.
associated :: Kneading -> Maybe (Kneading, Kneading)
associated (StarPeriodic k) = Just (Periodic a, Periodic abar)
  where
    n = genericLength k
    divisors = [ m | m <- [1 .. n], n `mod` m == 0 ]
    abar = head . filter (and . zipWith (==) a' . cycle) . map (`genericTake` a') $ divisors
    (a, a') = if ((n `elem`) . internalToList) `fmap` internalAddress (Periodic a1) == Just True then (a1, a2) else (a2, a1)
    a1 = map (\s -> case s of Star -> Zero ; t -> t) k
    a2 = map (\s -> case s of Star -> One  ; t -> t) k
associated _ = Nothing

-- | The upper associated kneading sequence.
upper :: Kneading -> Maybe Kneading
upper = fmap fst . associated

-- | The lower associated kneading sequence.
lower :: Kneading -> Maybe Kneading
lower = fmap snd . associated

-- | Angled internal addresses have angles between each integer in an
--   internal address.
data AngledInternalAddress
  = Unangled Integer
  | Angled Integer Angle AngledInternalAddress
  deriving (Read, Show, Eq, Ord, Data, Typeable)

-- | Angled internal address as a string.
prettyAngledInternalAddress :: AngledInternalAddress -> String
prettyAngledInternalAddress (Unangled n) = show n
prettyAngledInternalAddress (Angled n r a)
    | r /= 1/2  = show n ++ " " ++ show (numerator r) ++ "/" ++ show (denominator r) ++ " " ++ prettyAngledInternalAddress a
    | otherwise = show n ++ " " ++ prettyAngledInternalAddress a

-- | Builds a valid 'AngledInternalAddress' from a list, checking the
--   precondition that only the last 'Maybe Angle' should be 'Nothing',
--   and the 'Integer' must be strictly increasing.
angledFromList :: [(Integer, Maybe Angle)] -> Maybe AngledInternalAddress
angledFromList = fromList' 0
  where
    fromList' x [(n, Nothing)] | n > x = Just (Unangled n)
    fromList' x ((n, Just r) : xs) | n > x && 0 < r && r < 1 = Angled n r `fmap` fromList' n xs
    fromList' _ _ = Nothing

unsafeAngledFromList :: [(Integer, Maybe Angle)] -> AngledInternalAddress
unsafeAngledFromList = fromList' 0
  where
    fromList' x [(n, Nothing)] | n > x = Unangled n
    fromList' x ((n, Just r) : xs) | n > x && 0 < r && r < 1 = Angled n r (fromList' n xs)
    fromList' _ _ = error "Fractal.Mandelbrot.Address.unsafeAngledFromList"

-- | Convert an 'AngledInternalAddress' to a list.
angledToList :: AngledInternalAddress -> [(Integer, Maybe Angle)]
angledToList (Unangled n) = [(n, Nothing)]
angledToList (Angled n r a) = (n, Just r) : angledToList a

denominators :: InternalAddress -> Kneading -> [Integer]
denominators a v = denominators' (internalToList a)
  where
    denominators' (s0:ss@(s1:_)) =
      let rr = r s0 s1
      in  (((s1 - rr) `div` s0) + if s0 `elem` takeWhile (<= s0) (orbit p rr) then 1 else 2) : denominators' ss
    denominators' _ = []
    r s s' = case s' `mod` s of
      0 -> s
      t -> t
    p = rho v

numerators :: Angle -> InternalAddress -> [Integer] -> [Integer]
numerators r a qs = zipWith num (internalToList a) qs
  where
    num s q = genericLength . filter (<= r) . map (genericIndex rs) $ [0 .. q - 2]
      where
        rs = iterate (foldr (.) id . genericReplicate s $ double) r

-- | The angled internal address corresponding to an external angle.
angledInternalAddress :: Angle -> Maybe AngledInternalAddress
angledInternalAddress r0 = do
  let r = wrap r0
      k = kneading r
  i <- internalAddress k
  let d = denominators i k
      n = numerators r i d
  return . unsafeAngledFromList . zip (internalToList i) . (++ [Nothing]) . map Just . zipWith (%) n $ d

-- | Split an angled internal address at the last island.
splitAddress :: AngledInternalAddress -> (AngledInternalAddress, [Angle])
splitAddress a =
  let (ps0, rs0) = unzip $ angledToList a
      ps1 = reverse ps0
      rs1 = reverse (Nothing : init rs0)
      prs1 = zip ps1 rs1
      f ((p, Just r):qrs@((q, _):_)) acc
        | p == denominator r * q = f qrs (r : acc)
      f prs acc = g prs acc
      g prs acc =
        let (ps2, rs2) = unzip prs
            ps3 = reverse ps2
            rs3 = reverse (Nothing : init rs2)
            prs3 = zip ps3 rs3
            aa = unsafeAngledFromList prs3
        in  (aa, acc)
  in  f prs1 []

-- | The inverse of 'splitAddress'.
joinAddress :: AngledInternalAddress -> [Angle] -> AngledInternalAddress
joinAddress (Unangled p) [] = Unangled p
joinAddress (Unangled p) (r:rs) = Angled p r (joinAddress (Unangled $ p * denominator r) rs)
joinAddress (Angled p r a) rs = Angled p r (joinAddress a rs)

-- | The period of an angled internal address.
addressPeriod :: AngledInternalAddress -> Integer
addressPeriod (Unangled p) = p
addressPeriod (Angled _ _ a) = addressPeriod a

-- | Discard angle information from an internal address.
stripAngles :: AngledInternalAddress -> InternalAddress
stripAngles = InternalAddress . map fst . angledToList

-- | The pair of external angles whose rays land at the root of the
--   hyperbolic component described by the angled internal address.
externalAngles :: AngledInternalAddress -> Maybe (Rational, Rational)
externalAngles = externalAngles' 1 (0, 1)

externalAngles' :: Integer -> (Rational, Rational) -> AngledInternalAddress -> Maybe (Rational, Rational)
externalAngles' p0 lohi a0@(Unangled p)
  | p0 /= p = case wakees lohi p of
      [lh] -> externalAngles' p lh a0
      _ -> Nothing
  | otherwise = Just lohi
externalAngles' p0 lohi a0@(Angled p r a)
  | p0 /= p = case wakees lohi p of
      [lh] -> externalAngles' p lh a0
      _ -> Nothing
  | otherwise = do
{-
      let num = numerator r
          den = denominator r
          q = p * den
          ws = wakees lohi q
          nums = [ num' | num' <- [ 1.. den - 1 ], let r' = num' % den, denominator r' == den ]
          nws, nnums :: Integer
          nws = genericLength ws
          nnums = genericLength nums
      guard (nws == nnums)
      i <- genericElemIndex num nums
      lh <- safeGenericIndex ws (i :: Integer)
      externalAngles' q lh a
-}
      let num = numerator r
          den = denominator r
          ws = wakees (0, 1) den
          nums = [ num' | num' <- [ 1.. den - 1 ], let r' = num' % den, denominator r' == den ]
          nws, nnums :: Integer
          nws = genericLength ws
          nnums = genericLength nums
      guard (nws == nnums)
      i <- genericElemIndex num nums
      (l,h) <- safeGenericIndex ws (i :: Integer)
      externalAngles' (p * den) (if p > 1 then (tune l lohi, tune h lohi) else (l, h)) a
wakees :: (Rational, Rational) -> Integer -> [(Rational, Rational)]
wakees (lo, hi) q =
  let gaps (l, h) n
        | n == 0 = [(l, h)]
--        | h - l < 1 % (2 ^ n - 1) = [(l, h)]
        | n > 0 = let gs = gaps (l, h) (n - 1)
                      cs = candidates n gs
                  in  accumulate cs gs
        | otherwise = error "Fractal.Mandelbrot.Address.gaps !(n >= 0)"
      candidates n gs =
        let den = 2 ^ n - 1
        in  [ r
            | (l, h) <- gs
            , num <- [ ceiling (l * fromInteger den)
                      .. floor (h * fromInteger den) ]
            , let r = num % den
            , l < r, r < h
            , period (kneading r) == Just n
            ]
      accumulate [] ws = ws
      accumulate (l : h : lhs) ws =
        let (ls, ms@((ml, _):_)) = break (l `inside`) ws
            (_s, (_, rh):rs) = break (h `inside`) ms
        in  ls ++ [(ml, l)] ++ accumulate lhs ((h, rh) : rs)
      accumulate _ _ = error "Fractal.Mandelbrot.Address.gaps !even"
      inside x (l, h) = l < x && x < h
  in  chunk2 . candidates q . gaps (lo, hi) $ (q - 1)

chunk2 :: [t] -> [(t, t)]
chunk2 [] = []
chunk2 (x:y:zs) = (x, y) : chunk2 zs
chunk2 _ = error "Fractal.Mandelbrot.Address.chunk2 !even"

genericElemIndex :: (Eq a, Integral b) => a -> [a] -> Maybe b
genericElemIndex _ [] = Nothing
genericElemIndex e (f:fs)
  | e == f = Just 0
  | otherwise = (1 +) `fmap` genericElemIndex e fs

safeGenericIndex :: Integral b => [a] -> b -> Maybe a
safeGenericIndex [] _ = Nothing
safeGenericIndex (x:xs) i
  | i < 0 = Nothing
  | i > 0 = safeGenericIndex xs (i - 1)
  | otherwise = Just x

-- | Parse an angle.
parseAngle :: String -> Maybe Angle
parseAngle s = case runP pFraction () "" s of
  Left _ -> Nothing
  Right f -> Just (unFraction f)

-- | Parse a list of angles.
parseAngles :: String -> Maybe [Angle]
parseAngles s = case runP (many pFraction) () "" s of
  Left _ -> Nothing
  Right fs -> Just (map unFraction fs)

-- | Parse a kneading element.
parseKnead :: String -> Maybe Knead
parseKnead s = case runP pKnead () "" s of
  Left _ -> Nothing
  Right k -> Just k

-- | Parse a non-aperiodic kneading sequence.
parseKneading :: String -> Maybe Kneading
parseKneading s = case runP pKneading () "" s of
  Left _ -> Nothing
  Right ks -> Just ks

-- | Parse an internal address.
parseInternalAddress :: String -> Maybe InternalAddress
parseInternalAddress s = case runP (many pNumber) () "" s of
  Left _ -> Nothing
  Right ns -> internalFromList (map unNumber ns)

-- | Parse an angled internal address, accepting some unambiguous
--   abbreviations.
parseAngledInternalAddress :: String -> Maybe AngledInternalAddress
parseAngledInternalAddress s = case runP parser () "" s of
  Left _ -> Nothing
  Right a -> Just a

data Token = Number Integer | Fraction Integer Integer

unFraction :: Token -> Angle
unFraction (Fraction t b) = t % b
unFraction _ = error "Fractal.Mandelbrot.Address.unFraction"

unNumber :: Token -> Integer
unNumber (Number n) = n
unNumber _ = error "Fractal.Mandelbrot.Address.unNumber"


type Parse t = ParsecT String () Identity t

parser :: Parse AngledInternalAddress
parser = do
  ts <- pTokens
  accum 1 ts
  where
    accum p [] = return $ Unangled p
    accum _ [Number n] = return $ Unangled n
    accum _ (Number n : ts@(Number _ : _)) = do
      a <- accum n ts
      return $ Angled n (1%2) a
    accum _ (Number n : Fraction t b : ts) = do
      a <- accum (n * b) ts
      return $ Angled n (t%b) a
    accum p (Fraction t b : ts) = do
      a <- accum (p * b) ts
      return $ Angled p (t % b) a

pTokens :: Parse [Token]
pTokens = do
  _ <- pOptionalSpace
  ts <- pToken `sepEndBy` pSpace
  eof
  return ts

pToken :: Parse Token
pToken = choice [ try pFraction, pNumber ]

pFraction :: Parse Token
pFraction = do
  Number top <- pNumber
  _ <- pOptionalSpace
  _ <- string "/"
  _ <- pOptionalSpace
  Number bottom <- pNumber
  guard  $ top < bottom
  return $ Fraction top bottom

pNumber :: Parse Token
pNumber = do
  n <- foldl (\x y -> 10 * x + y) 0 `fmap` map (toInteger . digitToInt) `fmap` many1 digit
  guard  $ 0 < n
  return $ Number n

pSpace :: Parse [String]
pSpace = many1 (string " ")

pOptionalSpace :: Parse [String]
pOptionalSpace = many (string " ")

pKnead :: Parse Knead
pKnead = choice [ string "0" >> return Zero, string "1" >> return One, string "*" >> return Star ]

pKneading :: Parse Kneading
pKneading = do
  pre <- many  pKnead
  _ <- string "("
  per <- many1 pKnead
  _ <- string ")"
  return $ case (null pre, last per) of
    (False, _)   -> PrePeriodic pre per
    (True, Star) -> StarPeriodic per
    _            -> Periodic per
