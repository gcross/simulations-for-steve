-- @+leo-ver=4-thin
-- @+node:gcross.20100505233148.1260:@thin simulate.hs
-- @@language Haskell

-- @<< Language extensions >>
-- @+node:gcross.20091120111528.1238:<< Language extensions >>
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE UnicodeSyntax #-}
-- @-node:gcross.20091120111528.1238:<< Language extensions >>
-- @nl

-- @<< Import needed modules >>
-- @+node:gcross.20091120111528.1235:<< Import needed modules >>
import Acme.Dont

import Control.Monad
import Control.Monad.Error
import Control.Monad.Trans

import Data.Complex
import Data.ConfigFile
import Data.Int
import Data.IORef
import Data.UUID
import Data.Vec ((:.)(..))
import Data.Vec.Nat

import System.Environment
import System.Posix.Clock
import System.Exit

import Text.Printf

import VMPS.Algorithms
import VMPS.Database
import VMPS.EnergyMinimizationChain
import VMPS.Operators
import VMPS.Operators.Dimensions
import VMPS.Models
import VMPS.States
import VMPS.Tensors

import Debug.Trace
-- @-node:gcross.20091120111528.1235:<< Import needed modules >>
-- @nl

-- @+others
-- @+node:gcross.20100505233148.1270:Values
-- @+node:gcross.20100505233148.1271:i
i :: Complex Double
i = 0 :+ 1
-- @-node:gcross.20100505233148.1271:i
-- @-node:gcross.20100505233148.1270:Values
-- @+node:gcross.20091120111528.1234:Operator tensors
phi_block_length = 4
lambda_block_length = 9

phi_block_start = 2
lambda_block_start = phi_block_start + phi_block_length
final = lambda_block_start + lambda_block_length

-- @+others
-- @+node:gcross.20100505233148.1275:b - 3
-- @+others
-- @+node:gcross.20100505233148.1269:H_phi
b_H_phi :: Int → OperatorSiteSpecification N3
b_H_phi block_start =
    [-- k = 1
     (1 ⇨ block_start+0)
     .
     -- s = 3/4
     ( 3/4 *:)
     .
     -- A = diag([1 1 0])
     SingleSiteOperator $
        (1 :. 0 :. 0 :. ()) :.
        (0 :. 1 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. ()) :.
                            ()
    ,-- k = 2
     (1 ⇨ block_start+1)
     .
     -- s = 1/4
     ( 1/4 *:)
     .
     -- A = diag([1 -1 0])
     SingleSiteOperator $
        (1 :.   0  :. 0 :. ()) :.
        (0 :. (-1) :. 0 :. ()) :.
        (0 :.   0  :. 0 :. ()) :.
                              ()
    ,-- k = 3
     (1 ⇨ block_start+2)
     .
     -- s = 1/4
     ( 1/4 *:)
     .
     -- A = accumarray({2 1},1,[3 3]) + accumarray({1 2},1,[3 3])
     SingleSiteOperator $
        (0 :. 1 :. 0 :. ()) :.
        (1 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. ()) :.
                            ()
    ,-- k = 4
     (1 ⇨ block_start+3)
     .
     -- s = -1/4
     (-1/4 *:)
     .
     -- A = i*( accumarray({1 2},1,[3 3]) - accumarray({2 1},1,[3 3]) )
     SingleSiteOperator $
        (  0 :. i :. 0 :. ()) :.
        ((-i):. 0 :. 0 :. ()) :.
        (  0 :. 0 :. 0 :. ()) :.
                               ()
    ]
-- @nonl
-- @-node:gcross.20100505233148.1269:H_phi
-- @+node:gcross.20100506124738.1280:H_lambda
b_H_lambda :: Int → Int → OperatorSiteSpecification N3
b_H_lambda block_start final =
    [-- k = 1, B = diag([0 0 1])
     (block_start+0 ⇨ final) . SingleSiteOperator $
        (0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 1 :. ()) :.
                            ()
    ,-- k = 2, B = diag([1 0 0])
     (block_start+1 ⇨ final) . SingleSiteOperator $
        (1 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. ()) :.
                            ()
    ,-- k = 3, B = diag([0 1 0])
     (block_start+2 ⇨ final) . SingleSiteOperator $
        (0 :. 0 :. 0 :. ()) :.
        (0 :. 1 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. ()) :.
                            ()
    ,-- k = 4, B = accumarray({2 3},1,[3,3]) + accumarray({3 2},1,[3,3])
     (block_start+3 ⇨ final) . SingleSiteOperator $
        (0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 1 :. ()) :.
        (0 :. 1 :. 0 :. ()) :.
                            ()
    ,-- k = 5, B = i*( accumarray({2 3},1,[3,3]) - accumarray({3 2},1,[3,3]) )
     (block_start+4 ⇨ final) . SingleSiteOperator $
        (0 :.  0 :. 0 :. ()) :.
        (0 :.  0 :. i :. ()) :.
        (0 :.(-i):. 0 :. ()) :.
                             ()
    ,-- k = 6, B = accumarray({3 1},1,[3,3]) + accumarray({1 3},1,[3,3])
     (block_start+5 ⇨ final) . SingleSiteOperator $
        (0 :. 0 :. 1 :. ()) :.
        (0 :. 0 :. 0 :. ()) :.
        (1 :. 0 :. 0 :. ()) :.
                            ()
    ,-- k = 7, B = i*( accumarray({3 1},1,[3,3]) - accumarray({1 3},1,[3,3]) )
     (block_start+6 ⇨ final) . SingleSiteOperator $
        (0 :. 0 :.(-i):. ()) :.
        (0 :. 0 :.  0 :. ()) :.
        (i :. 0 :.  0 :. ()) :.
                             ()
    ,-- k = 8, B = accumarray({2 1},1,[3,3]) + accumarray({1 2},1,[3,3])
     (block_start+7 ⇨ final) . SingleSiteOperator $
        (0 :. 1 :. 0 :. ()) :.
        (1 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. ()) :.
                            ()
    ,-- k = 9, B = i*( accumarray({2 1},1,[3,3]) - accumarray({1 2},1,[3,3]) )
     (block_start+8 ⇨ final) . SingleSiteOperator $
        (0 :.(-i):. 0 :. ()) :.
        (i :.  0 :. 0 :. ()) :.
        (0 :.  0 :. 0 :. ()) :.
                            ()
    ]
-- @nonl
-- @-node:gcross.20100506124738.1280:H_lambda
-- @-others

b_operator_tensor = makeOperatorSiteTensorFromSpecification final final $
    [(1 ⇨ 1) identity, (final ⇨ final) identity]
    ++
    b_H_phi phi_block_start
    ++
    b_H_lambda lambda_block_start final
-- @nonl
-- @-node:gcross.20100505233148.1275:b - 3
-- @+node:gcross.20100505233148.1276:q - 5
-- @+others
-- @+node:gcross.20100505233148.1274:H_phi
q_H_phi :: Int → Int → OperatorSiteSpecification N5
q_H_phi block_start final =
    [-- k = 1, B = diag([1 1 0 0 0])
     (block_start+0 ⇨ final) . SingleSiteOperator $
        (1 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 1 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
                                      ()
    ,-- k = 2, B = diag([1 -1 0 0 0])
     (block_start+1 ⇨ final) . SingleSiteOperator $
        (1 :.  0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :.(-1):. 0 :. 0 :. 0 :. ()) :.
        (0 :.  0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :.  0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :.  0 :. 0 :. 0 :. 0 :. ()) :.
                                         ()
    ,-- k = 3, B = accumarray({1 2},1,[5 5]) + accumarray({2 1},1,[5 5])
     (block_start+2 ⇨ final) . SingleSiteOperator $
        (0 :. 1 :. 0 :. 0 :. 0 :. ()) :.
        (1 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
                                      ()
    ,-- k = 4, B = i*( accumarray({2 1},1,[5 5]) - accumarray({1 2},1,[5 5]) )
     (block_start+3 ⇨ final) . SingleSiteOperator $
        (0 :.(-i):. 0 :. 0 :. 0 :. ()) :.
        (i :.  0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :.  0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :.  0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :.  0 :. 0 :. 0 :. 0 :. ()) :.
                                      ()
    ]
-- @nonl
-- @-node:gcross.20100505233148.1274:H_phi
-- @+node:gcross.20100505233148.1277:H_lambda
q_H_lambda :: Complex Double → Int → OperatorSiteSpecification N5
q_H_lambda lambda block_start =
    [-- k = 1
     (1 ⇨ block_start+0)
     .
     times_common_factor
     .
     -- A = diag([0 0 0 0 1])
     SingleSiteOperator $
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 1 :. ()) :.
                                      ()
    ,-- k = 2
     (1 ⇨ block_start+1)
     .
     times_common_factor
     .
     -- s = L^2/2
     (lambda_squared/2 *:)
     .
     -- A = diag([0 0 0 1 0])
     SingleSiteOperator $
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 1 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
                                      ()
    ,-- k = 3
     (1 ⇨ block_start+2)
     .
     times_common_factor
     .
     -- s = L^2/2
     (lambda_squared/2 *:)
     .
     -- A = diag([0 0 1 0 0])
     SingleSiteOperator $
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 1 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
                                      ()
    ,-- k = 4
     (1 ⇨ block_start+3)
     .
     times_common_factor
     .
     -- s = -(L/sqrt(2))/2
     (-lambda/(2*sqrt 2) *:)
     .
     -- A = accumarray({3 5},1,[5,5]) + accumarray({5 3},1,[5,5])
     SingleSiteOperator $
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 1 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 1 :. 0 :. 0 :. ()) :.
                                      ()
    ,-- k = 5
     (1 ⇨ block_start+4)
     .
     times_common_factor
     .
     -- s = (L/sqrt(2))/2
     (lambda/(2*sqrt 2) *:)
     .
     -- A = i*( accumarray({3 5},1,[5,5]) - accumarray({5 3},1,[5,5]) )
     SingleSiteOperator $
        (0 :. 0 :.  0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :.  0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :.  0 :. 0 :. i :. ()) :.
        (0 :. 0 :.  0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :.(-i):. 0 :. 0 :. ()) :.
                                        ()
    ,-- k = 6
     (1 ⇨ block_start+5)
     .
     times_common_factor
     .
     -- s = (L/sqrt(2))/2
     (lambda/(2*sqrt 2) *:)
     .
     -- A = accumarray({4 5},1,[5,5]) + accumarray({5 4},1,[5,5])
     SingleSiteOperator $
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 1 :. ()) :.
        (0 :. 0 :. 0 :. 1 :. 0 :. ()) :.
                                      ()
    ,-- k = 7
     (1 ⇨ block_start+6)
     .
     times_common_factor
     .
     -- s = (L/sqrt(2))/2
     (lambda/(2*sqrt 2) *:)
     .
     -- i*( accumarray({4 5},1,[5,5]) - accumarray({5 4},1,[5,5]) )
     SingleSiteOperator $
        (0 :. 0 :. 0 :.  0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :.  0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :.  0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :.  0 :. i :. ()) :.
        (0 :. 0 :. 0 :.(-i):. 0 :. ()) :.
                                       ()
    ,-- k = 8
     (1 ⇨ block_start+7)
     .
     times_common_factor
     .
     -- s = -(L^2/2)/2
     (-lambda_squared/4 *:)
     .
     -- A = accumarray({3 4},1,[5,5]) + accumarray({4 3},1,[5,5])
     SingleSiteOperator $
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 1 :. 0 :. ()) :.
        (0 :. 0 :. 1 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
                                      ()
    ,-- k = 9
     (1 ⇨ block_start+8)
     .
     times_common_factor
     .
     -- s = (L^2/2)/2
     (lambda_squared/4 *:)
     .
     -- A = i*( accumarray({3 4},1,[5,5]) - accumarray({4 3},1,[5,5]) );
     SingleSiteOperator $
        (0 :. 0 :.  0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :.  0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :.  0 :. i :. 0 :. ()) :.
        (0 :. 0 :.(-i):. 0 :. 0 :. ()) :.
        (0 :. 0 :.  0 :. 0 :. 0 :. ()) :.
                                       ()
    ]
  where
    lambda_squared = lambda*lambda
    times_common_factor = (1/(1+lambda_squared) *:)
-- @nonl
-- @-node:gcross.20100505233148.1277:H_lambda
-- @+node:gcross.20100506160128.1281:H_in
q_H_in :: Int → OperatorSiteSpecification N5
q_H_in final =
    [-- Hin = diag( [0 1 0 0 0] );
     (1 ⇨ final) . SingleSiteOperator $
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 1 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
        (0 :. 0 :. 0 :. 0 :. 0 :. ()) :.
                                      ()
    ]
-- @nonl
-- @-node:gcross.20100506160128.1281:H_in
-- @+node:gcross.20100506160128.1287:H_U
q_H_U :: Int → OperatorSiteSpecification N5
q_H_U final =
    [-- Hin = kron( [1 -1;-1 1]/2 , eye(2) );  HU(dq,dq) = 0
     (1 ⇨ final) . SingleSiteOperator $
        (  1 :.  0 :.(-1):.  0 :. 0 :. ()) :.
        (  0 :.  1 :.  0 :.(-1):. 0 :. ()) :.
        ((-1):.  0 :.  1 :.  0 :. 0 :. ()) :.
        (  0 :.(-1):.  0 :.  1 :. 0 :. ()) :.
        (  0 :.  0 :.  0 :.  0 :. 0 :. ()) :.
                                           ()
    ]
-- @nonl
-- @-node:gcross.20100506160128.1287:H_U
-- @-others

first_operator_tensor lambda = makeOperatorSiteTensorFromSpecification 1 final $
    [(1 ⇨ 1) identity]
    ++
    q_H_U final
    ++
    q_H_in final
    ++
    q_H_lambda lambda lambda_block_start

q_operator_tensor lambda = makeOperatorSiteTensorFromSpecification final final $
    [(1 ⇨ 1) identity, (final ⇨ final) identity]
    ++
    q_H_U final
    ++
    q_H_phi phi_block_start final
    ++
    q_H_lambda lambda lambda_block_start
-- @nonl
-- @-node:gcross.20100505233148.1276:q - 5
-- @+node:gcross.20100506160128.1282:end - 2
last_operator_tensor = makeOperatorSiteTensorFromSpecification final 1 $    
    [(final ⇨ 1) identity]
    ++
    end_H_phi phi_block_start
-- @nonl
-- @+node:gcross.20100506160128.1285:H_phi
end_H_phi :: Int → OperatorSiteSpecification N2
end_H_phi block_start =
    [-- k = 1, B = diag([1 1])
     (block_start+0 ⇨ 1) . SingleSiteOperator $
        (1 :. 0 :. ()) :.
        (0 :. 1 :. ()) :.
                       ()
    ,-- k = 2, B = diag([1 -1)
     (block_start+1 ⇨ 1) . SingleSiteOperator $
        (1 :.  0 :. ()) :.
        (0 :.(-1):. ()) :.
                          ()
    ,-- k = 3, B = accumarray({1 2},1,[2 2]) + accumarray({2 1},1,[2 2])
     (block_start+2 ⇨ 1) . SingleSiteOperator $
        (0 :. 1 :. ()) :.
        (1 :. 0 :. ()) :.
                       ()
    ,-- k = 4, B = i*( accumarray({2 1},1,[2 2]) - accumarray({1 2},1,[2 2]) )
     (block_start+3 ⇨ 1) . SingleSiteOperator $
        (0 :.(-i):. ()) :.
        (i :.  0 :. ()) :.
                        ()
    ]
-- @nonl
-- @-node:gcross.20100506160128.1285:H_phi
-- @-node:gcross.20100506160128.1282:end - 2
-- @-others

makeModelOperatorSiteTensors lambda number_of_middle_qb_pairs
 | number_of_middle_qb_pairs < 0 =
    error "Negative number of middle q-b pairs specified!"
 | otherwise =
    [first_operator_tensor lambda,b_operator_tensor]
    ++
    (concat . replicate number_of_middle_qb_pairs)
        [q_operator_tensor lambda,b_operator_tensor]
    ++
    [last_operator_tensor]
-- @-node:gcross.20091120111528.1234:Operator tensors
-- @+node:gcross.20091120111528.1237:analyzeTrialEnergies
data TrialAnalysis = TrialDidBetter | TrialDidWorse | TrialDidTheSame

analyzeTrialEnergies tolerance best_energy trial_energy
    | best_energy - trial_energy > tolerance = TrialDidBetter
    | trial_energy - best_energy > tolerance = TrialDidWorse
    | otherwise = TrialDidTheSame
-- @-node:gcross.20091120111528.1237:analyzeTrialEnergies
-- @+node:gcross.20091120111528.1236:main
main = do
    args ← getArgs

    let lambda = (:+ 0) . read $ args !! 0
        number_of_middle_qb_pairs = read $ args !! 1
        operator_site_tensors = makeModelOperatorSiteTensors lambda number_of_middle_qb_pairs
        bandwidth_increment = 2
        initial_bandwidth = 2
        bandwidth_increase_energy_change_convergence_criterion = 1e-4
        multisweep_energy_change_convergence_criterion = 1e-4
        level_similarity_tolerance = 1e-3
        eigensolver_tolerance = 1e-10
        maximum_allowed_eigensolver_iterations = 1000

    putStrLn $
        printf "Running simulation for lambda = %s with %i sites..."
            (show . realPart $ lambda)
            (2*number_of_middle_qb_pairs+3)

    -- @    << Define callbacks >>
    -- @+node:gcross.20091205211300.1710:<< Define callbacks >>
    next_bandwidth_ref ← newIORef initial_bandwidth
    level_number_ref ← newIORef 1

    let getHeading = liftM (printf "LEVEL %i: ") (readIORef level_number_ref :: IO Int)
        callback_to_decide_whether_to_declare_victory_with_trial chain = do
            heading ← getHeading
            putStrLn $ heading ++ " energy = " ++ (show . chainEnergy $ chain)
            level_number ← readIORef level_number_ref
            let new_level_number = level_number + 1
            putStrLn $ printf "Now starting on level %i... (bandwidth=2 sweeps will not be displayed)" new_level_number
            writeIORef level_number_ref new_level_number
            alwaysDeclareVictory chain
        callback_to_increase_bandwidth chain = do
            next_bandwidth ← readIORef next_bandwidth_ref
            writeIORef next_bandwidth_ref (next_bandwidth+bandwidth_increment)
            increaseChainBandwidth next_bandwidth chain
        callback_after_each_sweep victory_flag latest_chain = do
            heading ← getHeading
            next_bandwidth ← readIORef next_bandwidth_ref
            let current_bandwidth = next_bandwidth-bandwidth_increment
            -- unless (current_bandwidth <= 2) $
            putStrLn $ heading ++ (printf "bandwidth = %i, sweep energy = %f" current_bandwidth (chainEnergy latest_chain) )
    -- @nonl
    -- @-node:gcross.20091205211300.1710:<< Define callbacks >>
    -- @nl

    -- @    << Run simulation >>
    -- @+node:gcross.20091202133456.1303:<< Run simulation >>
    (energies,_,_) ←
        fmap unzip3 $ solveForMultipleLevelsWithCallbacks
            callback_to_decide_whether_to_declare_victory_with_trial
            (newChainCreator
                (\number_of_projectors ->
                    let starting_bandwidth = initial_bandwidth+number_of_projectors
                    in writeIORef next_bandwidth_ref (starting_bandwidth+bandwidth_increment)
                       >>
                       return starting_bandwidth
                )
                operator_site_tensors
            )
            callback_to_increase_bandwidth
            callback_after_each_sweep
            ignoreSiteCallback
            bandwidth_increase_energy_change_convergence_criterion
            multisweep_energy_change_convergence_criterion
            eigensolver_tolerance
            maximum_allowed_eigensolver_iterations
            3
            []
    -- @-node:gcross.20091202133456.1303:<< Run simulation >>
    -- @nl

    putStrLn ""
    putStrLn "The energy levels are:"
    forM_ energies $ \energy → do
        putStr "\t"
        putStrLn . show $ energy

    let ground_energy:excited_energy:_ = energies
        energy_gap = excited_energy - ground_energy

    putStrLn ""
    putStrLn $ "The gap is " ++ show energy_gap

    TimeSpec time_in_seconds _ ← getTime ProcessCPUTime

    putStrLn $ "The elapsed CPU time for this run was " ++ show time_in_seconds ++ " seconds."
-- @nonl
-- @-node:gcross.20091120111528.1236:main
-- @-others
-- @-node:gcross.20100505233148.1260:@thin simulate.hs
-- @-leo
