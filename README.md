Heat_Treansfer_2D_finite_analysis

This repository contains a MATLAB implementation for solving 2D steady-state heat conduction problems using the Finite Difference Method (FDM). The solver uses Gauss-Seidel iteration to compute temperature distributions across an L-shaped domain, with support for linearized radiation boundary conditions alongside standard convective and fixed-temperature conditions.

Features
Finite difference discretization of the 2D steady-state heat equation over a non-rectangular (L-shaped) geometry
Gauss-Seidel iterative solver for the resulting linear system
Linearized treatment of radiative boundary conditions for improved convergence
Configurable boundary conditions (Dirichlet, Neumann/convective, radiative)
Temperature field visualization (contour/surface plots)
How It Works

The domain is discretized into a grid, and an energy balance is applied at each node to generate a system of algebraic equations. Boundary nodes account for the specified conditions (fixed temperature, convection, or radiation), with radiation terms linearized to keep the system solvable within the iterative Gauss-Seidel framework. The solver iterates until the temperature field converges within a specified tolerance.
