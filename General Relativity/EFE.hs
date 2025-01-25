-- Einsten Field Equations 

{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE AllowAmbiguousTypes #-}

module EinsteinFieldEquations where

import Data.Tensor.TypeLevel
import Data.Tensor.TypeLevel.TH
import Data.Coerce
import GHC.TypeLits
import Numeric.LinearAlgebra.Static (Matrix, V4, diag, inv, (><), det, vector, matrix, size)
import qualified Numeric.LinearAlgebra.Static as HMatrix
import Data.Maybe (fromMaybe)

type Dimension = 4

type SpacetimeMatrix = Matrix Dimension Dimension Double

type SpacetimeVector = V4 Double

type MetricTensor = SpacetimeMatrix

type InverseMetricTensor = SpacetimeMatrix

minkowskiMetric :: MetricTensor
minkowskiMetric = diag $ vector [-1, 1, 1, 1]

inverseMetric :: MetricTensor -> InverseMetricTensor
inverseMetric g = fromMaybe (error "Singular metric tensor") $ inv g

christoffelSymbols :: MetricTensor -> [[[Double]]]
christoffelSymbols g =
  let gInv = inverseMetric g
      dim = 4
      d_sigma_g_munu :: Int -> Int -> Int -> Double
      d_sigma_g_munu _ _ _ = 0

  in  [[[ christoffelSymbolComponent gInv d_sigma_g_munu lambda mu nu
         | nu <- [0..dim-1] ]
       | mu <- [0..dim-1] ]
     | lambda <- [0..dim-1] ]

  where
    christoffelSymbolComponent :: InverseMetricTensor -> (Int -> Int -> Int -> Double) -> Int -> Int -> Int -> Double
    christoffelSymbolComponent gInv' d_sigma_g_munu' lambda' mu' nu' =
      let sumTerm = sum [ gInv' ! (lambda', sigma') * ( term1 + term2 - term3 )
                        | sigma' <- [0..3]
                        , let term1 = 0.5 * d_sigma_g_munu' sigma' mu' nu'
                              term2 = 0.5 * d_sigma_g_munu' mu' nu' sigma'
                              term3 = 0.5 * d_sigma_g_munu' nu' sigma' mu'
                        ]
      in  sumTerm


riemannTensor :: MetricTensor -> [[[[Double]]]]
riemannTensor g =
  let christoffel = christoffelSymbols g
      dim = 4
      d_mu_christoffel_nu_rho_sigma :: Int -> Int -> Int -> Int -> Double
      d_mu_christoffel_nu_rho_sigma _ _ _ _ = 0

  in  [[[[ riemannTensorComponent christoffel d_mu_christoffel_nu_rho_sigma lambda sigma mu nu
          | nu <- [0..dim-1] ]
        | mu <- [0..dim-1] ]
      | sigma <- [0..dim-1] ]
    | lambda <- [0..dim-1] ]

  where
    riemannTensorComponent :: [[[Double]]] -> (Int -> Int -> Int -> Int -> Double) -> Int -> Int -> Int -> Int -> Double
    riemannTensorComponent christoffel' d_mu_christoffel_nu_rho_sigma' lambda' sigma' mu' nu' =
      let term1 = d_mu_christoffel_nu_rho_sigma' mu' lambda' sigma' nu'
          term2 = d_mu_christoffel_nu_rho_sigma' nu' lambda' sigma' mu'
          term3Sum = sum [ christoffel' !! lambda' !! rho' !! mu' * christoffel' !! rho' !! sigma' !! nu' | rho' <- [0..3] ]
          term4Sum = sum [ christoffel' !! lambda' !! rho' !! nu' * christoffel' !! rho' !! sigma' !! mu' | rho' <- [0..3] ]
      in  term1 - term2 + term3Sum - term4Sum


ricciTensor :: MetricTensor -> MetricTensor
ricciTensor g =
  let riemann = riemannTensor g
      dim = 4
  in  matrix dim dim $ \ (mu, nu) ->
        sum [ riemann !! rho !! mu !! rho !! nu | rho <- [0..dim-1] ]


ricciScalar :: MetricTensor -> Double
ricciScalar g =
  let ricci = ricciTensor g
      gInv = inverseMetric g
      dim = 4
  in  sum [ gInv ! (mu, nu) * (ricci ! (mu, nu)) | mu <- [0..dim-1], nu <- [0..dim-1] ]


einsteinTensor :: MetricTensor -> MetricTensor
einsteinTensor g =
  let ricci = ricciTensor g
      ricciS = ricciScalar g
      dim = 4
  in  matrix dim dim $ \ (mu, nu) ->
        ricci ! (mu, nu) - 0.5 * ricciS * (g ! (mu, nu))


vacuumEnergyMomentumTensor :: MetricTensor -> MetricTensor
vacuumEnergyMomentumTensor _ = HMatrix.zero


checkEinsteinFieldEquations :: MetricTensor -> MetricTensor -> Double -> Bool
checkEinsteinFieldEquations g energyMomentumTensor tolerance =
  let einsteinT = einsteinTensor g
      constantFactor = 8 * pi
      rhs = constantFactor * energyMomentumTensor g
      lhs = einsteinT

      dimension = 4
      areComponentsClose :: Bool
      areComponentsClose = and [ abs (lhs ! (mu, nu) - rhs ! (mu, nu)) < tolerance
                                | mu <- [0..dimension-1], nu <- [0..dimension-1] ]

  in  areComponentsClose


main :: IO ()
main = do
  putStrLn "Einstein Field Equations in Haskell (Simplified Illustration)"
  putStrLn "\nUsing Minkowski Metric (flat spacetime):"
  let gMinkowski = minkowskiMetric

  putStrLn "\nCalculating Christoffel Symbols (for Minkowski):"
  let christoffelMinkowski = christoffelSymbols gMinkowski
  putStrLn $ "Christoffel Symbols (first component Γ^0_00): " ++ show (christoffelMinkowski !! 0 !! 0 !! 0)

  putStrLn "\nCalculating Riemann Tensor (for Minkowski):"
  let riemannMinkowski = riemannTensor gMinkowski
  putStrLn $ "Riemann Tensor (first component R^0_000): " ++ show (riemannMinkowski !! 0 !! 0 !! 0 !! 0)

  putStrLn "\nCalculating Ricci Tensor (for Minkowski):"
  let ricciMinkowski = ricciTensor gMinkowski
  putStrLn $ "Ricci Tensor: " ++ show ricciMinkowski

  putStrLn "\nCalculating Ricci Scalar (for Minkowski):"
  let ricciScalarMinkowski = ricciScalar gMinkowski
  putStrLn $ "Ricci Scalar: " ++ show ricciScalarMinkowski

  putStrLn "\nCalculating Einstein Tensor (for Minkowski):"
  let einsteinMinkowski = einsteinTensor gMinkowski
  putStrLn $ "Einstein Tensor: " ++ show einsteinMinkowski

  putStrLn "\nChecking Einstein Field Equations in vacuum (T_μν = 0) for Minkowski:"
  let isEFEsatisfied = checkEinsteinFieldEquations gMinkowski vacuumEnergyMomentumTensor 1e-9
  putStrLn $ "Are EFE satisfied (within tolerance)? " ++ show isEFEsatisfied