# Numerical Methods

## Overview

This section covers the numerical methods and algorithms used to solve the discretized governing equations in CFD.

## Solution of Linear Systems

After discretization, the governing equations are typically reduced to a system of linear algebraic equations:

```
Ax = b
```

Where:
- A = coefficient matrix
- x = solution vector (unknowns)
- b = right-hand side vector

## Direct Methods

Direct methods solve the linear system exactly (within machine precision).

### Gaussian Elimination
- Systematic elimination of variables
- Computational cost: O(n³)
- Memory requirements: O(n²)
- Suitable for small systems only

### LU Decomposition
```
A = LU
```
- Decomposes A into lower (L) and upper (U) triangular matrices
- Efficient for multiple right-hand sides
- Cost: O(n³) for decomposition, O(n²) for back substitution

### Thomas Algorithm (TDMA)
- Special algorithm for tridiagonal matrices
- Very efficient: O(n)
- Common in 1D problems and line-by-line solvers

### Advantages of Direct Methods
- Exact solution (within roundoff error)
- Predictable computational cost
- No convergence issues

### Disadvantages of Direct Methods
- High memory requirements for large systems
- Expensive for sparse matrices
- Not practical for 3D problems (millions of unknowns)

## Iterative Methods

Iterative methods generate a sequence of approximations that converge to the solution.

### Point Iterative Methods

#### Jacobi Method
```
x_i^{k+1} = (b_i - Σ_{j≠i} a_{ij}x_j^k) / a_{ii}
```

- Simple to implement
- Slow convergence
- Can be parallelized easily

#### Gauss-Seidel Method
```
x_i^{k+1} = (b_i - Σ_{j<i} a_{ij}x_j^{k+1} - Σ_{j>i} a_{ij}x_j^k) / a_{ii}
```

- Uses most recent values
- Faster convergence than Jacobi
- Cannot be easily parallelized

#### Successive Over-Relaxation (SOR)
```
x_i^{k+1} = x_i^k + ω(x_i^{GS} - x_i^k)
```

- Relaxation factor: 1 < ω < 2
- Accelerates Gauss-Seidel
- Optimal ω problem-dependent

### Krylov Subspace Methods

More sophisticated iterative methods based on Krylov subspaces.

#### Conjugate Gradient (CG)
- For symmetric positive definite matrices
- Optimal in Krylov subspace
- Theoretically converges in n iterations
- Practical convergence much faster

#### Bi-Conjugate Gradient Stabilized (BiCGSTAB)
- For non-symmetric matrices
- More stable than standard BiCG
- Commonly used in CFD

#### Generalized Minimal Residual (GMRES)
- Minimizes residual in Krylov subspace
- Works for general matrices
- Memory requirements increase with iterations
- Often restarted after m iterations (GMRES(m))

### Multigrid Methods

Multigrid methods accelerate convergence by solving on multiple grid levels.

#### Two-Grid Cycle
1. **Relaxation**: Smooth error on fine grid
2. **Restriction**: Transfer residual to coarse grid
3. **Coarse grid correction**: Solve on coarse grid
4. **Prolongation**: Interpolate correction to fine grid
5. **Relaxation**: Post-smoothing on fine grid

#### V-Cycle and W-Cycle
- V-cycle: Go down to coarsest level, then up
- W-cycle: Multiple visits to intermediate levels
- Full Multigrid (FMG): Start from coarsest level

#### Advantages
- Optimal complexity: O(n)
- Very fast convergence
- Independent of grid size

#### Disadvantages
- Complex to implement
- Requires hierarchy of grids
- May require problem-specific tuning

## Pressure-Velocity Coupling

For incompressible flows, special algorithms are needed to couple pressure and velocity.

### SIMPLE Algorithm

**Semi-Implicit Method for Pressure-Linked Equations**

1. Guess pressure field p*
2. Solve momentum equations for velocity u*
3. Solve pressure correction equation
4. Correct pressure and velocity
5. Repeat until convergence

### SIMPLEC
- SIMPLE-Consistent
- Improved pressure correction
- Faster convergence than SIMPLE

### SIMPLER
- SIMPLE-Revised
- Solves pressure equation directly
- More robust for some problems

### PISO
- Pressure Implicit with Splitting of Operators
- Two or more corrector steps
- Better for transient problems
- Used in OpenFOAM

### Fractional Step Method
1. Solve momentum without pressure gradient (predictor)
2. Solve pressure Poisson equation (corrector)
3. Update velocity to satisfy continuity

## Convection-Diffusion Schemes

Schemes for discretizing convective terms.

### Upwind Schemes

#### First-Order Upwind
```
φ_face = φ_upwind
```
- Very stable
- Highly diffusive
- First-order accurate

#### Second-Order Upwind
```
φ_face = φ_upwind + ∇φ·Δx
```
- Less diffusive
- Second-order accurate
- May oscillate for sharp gradients

### Central Differencing
```
φ_face = (φ_i + φ_{i+1})/2
```
- Second-order accurate
- Unbounded (can oscillate)
- Good for diffusion-dominated flows

### Higher-Order Schemes

#### QUICK (Quadratic Upstream Interpolation)
- Third-order accurate
- Uses quadratic interpolation
- Better accuracy than second-order upwind

#### TVD Schemes (Total Variation Diminishing)
- Prevents spurious oscillations
- Preserves monotonicity
- Examples: MUSCL, van Leer, Superbee

### Blended Schemes
```
φ_face = γφ_upwind + (1-γ)φ_central
```
- Combine stability and accuracy
- γ = 0: Central (accurate but may oscillate)
- γ = 1: Upwind (stable but diffusive)

## Convergence Criteria

### Residual Monitoring
```
Residual = |Ax^k - b|
```

Typical convergence criteria:
- Relative residual < 10⁻³ to 10⁻⁶
- Monitor all equations (continuity, momentum, energy)

### Variable Monitoring
- Monitor solution variables (velocity, pressure, temperature)
- Check for steady-state (∂φ/∂t → 0)
- Verify global conservation (mass, energy)

### Under-Relaxation

Used to improve stability of iterative process:
```
φ^{k+1} = α·φ^{new} + (1-α)·φ^k
```

Where 0 < α ≤ 1:
- α = 1: No relaxation
- α < 1: Under-relaxation (more stable, slower convergence)

Typical values:
- Pressure: 0.3-0.7
- Momentum: 0.7-0.9
- Turbulence: 0.5-0.8

## Time Integration Methods

### Explicit Methods

#### Forward Euler
```
φ^{n+1} = φ^n + Δt·F(φ^n, t^n)
```
- First-order accurate
- Conditionally stable (CFL restriction)

#### Runge-Kutta Methods
- Higher-order accuracy
- Multiple stages per time step
- Popular: RK4 (fourth-order)

### Implicit Methods

#### Backward Euler
```
φ^{n+1} = φ^n + Δt·F(φ^{n+1}, t^{n+1})
```
- First-order accurate
- Unconditionally stable
- Requires solving non-linear system

#### Crank-Nicolson
```
φ^{n+1} = φ^n + Δt·[F(φ^{n+1}, t^{n+1}) + F(φ^n, t^n)]/2
```
- Second-order accurate
- Unconditionally stable
- Good balance of accuracy and stability

## Convergence Acceleration

### Residual Smoothing
- Reduces high-frequency errors
- Allows larger time steps
- Common in explicit schemes

### Local Time Stepping
- Different Δt for each cell
- Based on local CFL condition
- Accelerates convergence to steady-state

### Preconditioning
- Modify system to improve conditioning
- Reduces iteration count
- Examples: ILU, Jacobi, multigrid
