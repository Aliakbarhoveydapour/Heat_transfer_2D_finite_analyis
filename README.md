# 2D Steady-State Heat Conduction — Finite Difference / Gauss-Seidel Solver

A MATLAB solver for steady-state heat conduction in a 2D L-shaped composite
solid, discretized with the finite difference method (FDM) and solved
iteratively with Gauss-Seidel. The model combines five different boundary
condition types on a single 41-node mesh, including a **linearized radiation**
term at the hot face — a common simplification used to keep a radiation
boundary condition solvable within a linear iterative scheme.

## Problem Description

The geometry is an L-shaped 2D solid (e.g. a bracket/fin cross-section)
discretized into a 41-node square grid (Δx = Δy = 10 mm), with the following
boundary conditions:

| Boundary | Type | Condition |
|---|---|---|
| Bottom face | Dirichlet | Fixed temperature, T = 60 °C |
| Left face | Mixed | Convection (h₁ = 50 W/m²·K, T∞ = 200 °C) + radiation (ε = 0.8, T_sur = 150 °C) |
| Right face (thin base) | Convection | h₂ = 20 W/m²·K, T∞ = 25 °C |
| Top of thin base | Insulated | q″ = 0 |
| All other outer faces | Convection | h₃ = 25 W/m²·K, T∞ = 25 °C |

Material: k = 180 W/(m·K) (aluminum-like conductivity).

## Method

- **Discretization:** standard 2D FDM energy balances, with corner/edge nodes
  derived from the appropriate control-volume balance for each boundary type.
- **Radiation linearization:** the nonlinear radiation term
  q_rad = εσ(T⁴ - T_sur⁴) is rewritten as q_rad = h_r(T)(T_sur - T), with
  h_r(T) = εσ(T_sur² + T²)(T_sur + T) recomputed from the current node
  temperature each iteration. This keeps the system linear per Gauss-Seidel
  sweep while still capturing T⁴ physics.
- **Solver:** point Gauss-Seidel iteration, convergence tolerance 1e-6 K,
  max 50,000 iterations.
- **Output:** converged nodal temperatures (°C) printed to the console.

## Requirements

- MATLAB R2023b (no toolboxes required — plain scalar MATLAB)

## Usage

```matlab
% From the MATLAB command window / editor:
heat_transfer_fdm_gauss_seidel.m
```

The script runs standalone — no input files or function arguments needed.
Edit the parameters at the top of the script (`k`, `dx`, boundary
temperatures, film coefficients, emissivity) to adapt it to a different
geometry or material.

## Sample Output

```
Running Gauss-Seidel with linearised radiation
Converged in <N> iterations  (residual = <e> K)

Node     Temp (°C)
----     ---------
Node  1:   60.000 °C
Node  2:   60.000 °C
...
```

## Notes

- Node numbering follows a hand-drawn mesh map of the L-shaped domain
  (see comments in the script for how each node connects to its neighbors).
- This was developed as a heat transfer coursework project; it's shared here
  as a worked example of coupling FDM boundary-condition derivation with an
  iterative linear solver, including a nonlinear (radiative) BC handled via
  linearization.

## Author

Aliakbar Hoveydapour
