import Data.Ratio
import Control.Monad

vdc :: Integer -> Integer -> Rational 
vdc base = go 0 1
  where
    go vdc _ 0 = vdc
    go vdc denom n =
      let
        denom' = denom * base
        (n', r) = n `quotRem` base
      in go (vdc + r % denom') denom' n' 

main :: IO ()
main = forM_ [3,5] $ \base -> do
         print base 
         forM_ [0..10000] (print . vdc base)
