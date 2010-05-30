-- @+leo-ver=4-thin
-- @+node:gcross.20100530102143.1307:@thin plot-gap-versus-number-of-sites.hs
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
        lambda_list :: [String] <- liftIO getArgs
        forM lambda_list $ \lambda ->
            let sql_statement = "select number_of_sites, energy_gap from simulations where lambda = " ++ lambda ++ " order by number_of_sites asc;"
            in do
                liftIO . hPutStrLn stderr $ "> " ++ sql_statement
                (data_points :: [(Float,Float)]) <- doQuery (sql $ sql_statement) fetch2 []
                return (lambda,data_points)
    )
    >>=
    return . concatMap (\(color,(lambda,data_points)) ->
        [ (PlotStyle Lines (CustomStyle [LineTitle ("lambda = " ++ lambda), LineType color, PointType 2]),data_points)
        ]
    ) . zip [1..]
    >>=
    plotPathsStyle
        [Custom "logscale" ["xy"]
        ,Key (Just ["right","top"])
        ,Title "Energy Gap vs. Number of Sites"
        ,XLabel "Number of Sites"
        ,YLabel "Energy Gap"
        ]
-- @-node:gcross.20100530102143.1307:@thin plot-gap-versus-number-of-sites.hs
-- @-leo
