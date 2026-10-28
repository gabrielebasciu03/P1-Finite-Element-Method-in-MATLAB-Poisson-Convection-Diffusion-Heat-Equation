# P1-Finite-Element-Method-in-MATLAB-Poisson-Convection-Diffusion-Heat-Equation
From-scratch P1 Galerkin finite element solvers in MATLAB: Poisson with PCG, convection-diffusion with SUD stabilisation, and the heat equation with the θ-method.

## Overview

This repository contains three MATLAB projects for the Numerical Methods for Partial Differential Equations course (University of Padua, A.Y. 2025-2026). All of them build on a from-scratch P1 Galerkin finite element solver on triangular meshes (mesh0 to mesh4, with h halved at each level). The shared infrastructure covers mesh input, local basis coefficients, element-wise sparse assembly of the stiffness, mass, transport and stabilisation matrices, and Dirichlet conditions through the penalty method.

**Project 1: Poisson equation.** Solves −Δu = f with mixed Dirichlet and Neumann conditions on (−1,1)² and compares the result with a known exact solution. The linear system is solved with backslash and with PCG preconditioned by IC(0). The PCG iteration count roughly doubles at each refinement (8 to 109). The L² error converges with a global order of 2.10, matching the theoretical O(h²).

**Project 2: Convection-diffusion.** Solves −DΔu + β·∇u = 0 with β = (1,3) on the unit square, with a boundary layer near y = 1. Without stabilisation, the Galerkin solution oscillates when convection dominates (max u = 4.11 at D = 10⁻³). A scan over D finds the critical diffusion D ≈ 7·10⁻² and a grid Péclet number Pe_h ≈ 2. Streamline-Upwind Diffusion (SUD) stabilisation is then tested for several values of δ, showing the trade-off between suppressing oscillations and adding numerical diffusion.

**Project 3: Heat equation.** Solves ∂u/∂t = Δu + f with the θ-method (explicit Euler, Crank-Nicolson, implicit Euler) and a manufactured solution. The spatial order is about 1.95. Direct temporal convergence on a fixed mesh is masked by the spatial error floor, so self-convergence (Δt against Δt/2) is used, recovering orders of 1.02 for implicit Euler and 1.92 for Crank-Nicolson. For explicit Euler, the stability limit Δt ≤ 2/λmax(P⁻¹A) is computed from a generalised eigenvalue problem, and the number of time steps grows from 74 to 19,595 across the meshes.

**Skills demonstrated:** finite element method, sparse matrix assembly, iterative solvers and preconditioning, numerical stabilisation, time integration and stability analysis, convergence studies, MATLAB.

## Repository contents

| File | Description |
|---|---|
| `LAB1.m` | Project 1: Poisson solver, error analysis and PCG study |
| `LAB2.m` | Project 2: convection-diffusion solver with SUD stabilisation |
| `LAB3.m` | Project 3: heat equation with the θ-method and convergence studies |
| Helper functions | `input_data*.m`, `localBasis*.m`, `stiffBuild*.m`, `massBuild3.m`, `transportBuild.m`, `stabBuild2.m`, `rhsBuild*.m`, `imposeBC*.m`, `pcgSolve.m`, `solveParabolic.m` and related routines |
| `mesh0` to `mesh4` | Triangulations provided by the course (`triang.dat`, `xy.dat`, `dirnod.dat`, `neunod.dat`) |
| `Report.pdf` | Written report with figures and comments |

(Rename the files to match your repository.)

## Project 1: Poisson equation

Domain (−1,1)², exact solution u = sin(πx)cos(πy/2), homogeneous Dirichlet on y = ±1 and Neumann on x = ±1.

| Mesh | h | L² error | Ratio |
|---|---|---|---|
| mesh0 | 7.1113e-01 | 2.5124e-01 | – |
| mesh1 | 3.5557e-01 | 5.1095e-02 | 4.9171 |
| mesh2 | 1.7778e-01 | 1.1994e-02 | 4.2601 |
| mesh3 | 8.8892e-02 | 2.9470e-03 | 4.0699 |
| mesh4 | 4.4446e-02 | 7.3346e-04 | 4.0179 |

Global convergence order: 2.0956.

| Mesh | ‖uh − uh,PCG‖ | PCG iterations |
|---|---|---|
| mesh0 | 5.522e-09 | 8 |
| mesh1 | 1.235e-08 | 15 |
| mesh2 | 6.658e-09 | 29 |
| mesh3 | 7.295e-08 | 55 |
| mesh4 | 8.935e-08 | 109 |

## Project 2: Convection-diffusion

Unit square, β = (1,3), u = 1 on Γ1 and u = 0 on Γ2, computed on mesh2.

| Quantity | Value |
|---|---|
| Critical diffusion D_crit | 7.0170e-02 |
| Mesh size h | 8.8892e-02 |
| Pseudo-Péclet number Pe_h = \|β\|h / (2·D_crit) | 2.0030 |

| D | δ (SUD) | max(u) | min(u) |
|---|---|---|---|
| 1 | 0 (none) | 1.0000 | 0.0000 |
| 1e-1 | 0 (none) | 1.0000 | 0.0000 |
| 1e-3 | 0 (none) | 4.1053 | −1.2532 |
| 1e-3 | 5e-2 | 1.0226 | −0.0953 |
| 1e-3 | 1e-2 | 1.0526 | −0.1231 |
| 1e-3 | 1e-3 | 1.9670 | −0.3312 |

## Project 3: Heat equation

Domain (−1,1)², T = 1, manufactured solution u = e⁻ᵗ sin(πx)cos(πy/2).

Spatial convergence (Δt = 1e-4):

| Mesh | h | Error (θ=1) | Error (θ=0.5) |
|---|---|---|---|
| mesh0 | 5.0285e-01 | 1.0144e-01 | 1.0154e-01 |
| mesh1 | 2.5142e-01 | 3.0472e-02 | 3.0522e-02 |
| mesh2 | 1.2571e-01 | 7.9624e-03 | 7.9656e-03 |
| mesh3 | 6.2856e-02 | 2.0009e-03 | 2.0036e-03 |
| mesh4 | 3.1428e-02 | 4.9912e-04 | 5.0166e-04 |

Orders: 1.9527 (θ = 1) and 1.9511 (θ = 0.5).

Temporal self-convergence on mesh4:

| Δt | ‖u^Δt − u^(Δt/2)‖ (θ=1) | ‖u^Δt − u^(Δt/2)‖ (θ=0.5) |
|---|---|---|
| 0.200000 | 1.8007e-03 | 8.9911e-05 |
| 0.100000 | 8.5619e-04 | 2.2077e-05 |
| 0.050000 | 4.1647e-04 | 6.7938e-06 |
| 0.025000 | 2.0533e-04 | 3.1794e-06 |
| 0.012500 | 1.0194e-04 | 1.3407e-06 |
| 0.006250 | 5.0790e-05 | 1.2033e-07 |
| 0.003125 | 2.5350e-05 | 1.9774e-08 |

Orders: 1.0226 (θ = 1) and 1.9226 (θ = 0.5). The direct temporal orders (0.87 and 0.09) are distorted by the spatial error floor of mesh4.

Explicit Euler (θ = 0) at the stability limit:

| Mesh | h | Δt_stab | Steps | Error |
|---|---|---|---|---|
| mesh0 | 5.0285e-01 | 1.364e-02 | 74 | 1.1181e-01 |
| mesh1 | 2.5142e-01 | 3.427e-03 | 292 | 3.1230e-02 |
| mesh2 | 1.2571e-01 | 8.402e-04 | 1191 | 7.9923e-03 |
| mesh3 | 6.2856e-02 | 2.076e-04 | 4818 | 2.0091e-03 |
| mesh4 | 3.1428e-02 | 5.104e-05 | 19595 | 5.0296e-04 |

Order at the stability limit: 1.9551.

## Getting started

1. Place the mesh folders (`mesh0` to `mesh4`) and the helper functions in the MATLAB path.
2. Run `LAB1.m`, `LAB2.m` or `LAB3.m` from the MATLAB command window.
3. Each script prints the tables above and produces the convergence and solution plots.

`LAB3.m` uses `eigs`, which needs no extra toolbox, and `ichol` in Project 1 requires a recent MATLAB release.

## Limitations

- Project 2 has no exact solution, so it is validated through the maximum principle and qualitative behaviour only.
- The explicit Euler runs on the finest meshes need tens of thousands of steps, which is why stable schemes are preferred in practice.
- The penalty method with a very small ε can give misleading residuals with iterative solvers, so Project 2 uses ε = 10⁻¹⁰.

## Author

**Gabriele Basciu**, University of Padua, A.Y. 2025-2026
