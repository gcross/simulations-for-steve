-- @+leo-ver=4-thin
-- @+node:gcross.20100530003420.1952:@thin plot-gap-versus-lambda.hs
-- @@language Haskell

{-# LANGUAGE ScopedTypeVariables #-}

import Control.Arrow
import Control.Monad
import Control.Monad.Trans
import Database.Enumerator
import Database.PostgreSQL.Enumerator
import Data.Char
import Graphics.Gnuplot.Simple
import System
import System.IO
import System.Console.GetOpt
import Text.Printf

import VMPS.Database

main =
    makeConnection "reader"
    >>=
    flip withSession (do
        number_of_sites_list :: [Int] <- fmap (map read) $ liftIO getArgs
        forM number_of_sites_list $ \number_of_sites ->
            let sql_statement = "select lambda, energy_gap from simulations where number_of_sites = " ++ show number_of_sites ++ " order by lambda asc;"
            in do
                liftIO . hPutStrLn stderr $ "> " ++ sql_statement
                (data_points :: [(Float,Float)]) <- doQuery (sql $ sql_statement) fetch2 []
                return (number_of_sites,data_points)
    )
    >>=
    return . concatMap (\(color,(number_of_sites,data_points)) ->
        [ (PlotStyle Lines (CustomStyle [LineTitle ("number of sites = " ++ show number_of_sites), LineType color, PointType 2]),data_points)
        ]
    ) . zip [1..]
    >>=
    plotPathsStyle
        -- [Custom "logscale" ["xy"]
        [Key (Just ["right","top"])
        ,Title "Energy Gap vs. Lambda"
        ,XLabel "Lambda"
        ,YLabel "Energy Gap"
        ]
-- @-node:gcross.20100530003420.1952:@thin plot-gap-versus-lambda.hs
-- @-leo
