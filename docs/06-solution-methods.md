# Solution Methods

## Overview

This section covers the overall solution strategies and algorithms used in CFD to solve the complete system of governing equations.

## Segregated vs. Coupled Solvers

### Segregated Approach

Solves the governing equations sequentially:

1. Solve x-momentum equation for u
2. Solve y-momentum equation for v
3. Solve z-momentum equation for w
4. Solve pressure correction equation
5. Update velocities and pressure
6. Solve energy equation (if needed)
7. Solve turbulence equations (if needed)
8. Repeat until convergence

**Advantages:**
- Lower memory requirements
- More robust for certain problems
- Easier to implement
- Good for incompressible flows

**Disadvantages:**
- Slower convergence
- Strong coupling between equations
- Requires under-relaxation
- Many outer iterations needed

### Coupled Approach

Solves momentum and pressure equations simultaneously:

1. Solve momentum and continuity equations together
2. Solve energy equation (if needed)
3. Solve turbulence equations (if needed)
4. Repeat until convergence

**Advantages:**
- Faster convergence
- Better for compressible flows
- Fewer iterations
- Better for highly coupled problems

**Disadvantages:**
- Higher memory requirements
- More complex implementation
- Can be less robust
- Requires more sophisticated preconditioners

## Pressure-Based vs. Density-Based Solvers

### Pressure-Based Solvers

Use pressure as primary variable for incompressible/low-Mach flows.

**Algorithm:**
1. Solve momentum equations with guessed pressure
2. Solve pressure correction equation to satisfy continuity
3. Correct pressure and velocity
4. Update other variables

**Applications:**
- Incompressible flows (Ma < 0.3)
- Low-speed flows
- Natural convection
- Internal flows

**Examples:**
- SIMPLE family algorithms
- Fractional step methods
- Projection methods

### Density-Based Solvers

Use continuity equation with density as primary variable.

**Algorithm:**
1. Solve continuity, momentum, and energy equations coupled
2. Update density from equation of state
3. Update pressure from equation of state
4. Repeat until convergence

**Applications:**
- High-speed compressible flows
- Supersonic/hypersonic flows
- Aerodynamics with shock waves
- Gas dynamics

**Examples:**
- Jameson's scheme
- Roe's scheme
- AUSM family

## Pressure-Velocity Coupling Algorithms

### SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)

**Step-by-step procedure:**

1. **Set initial conditions**
   - Guess pressure field p*
   - Guess velocity field u*

2. **Solve momentum equations**
   ```
   a_P u*_P = Σ a_nb u*_nb + b + (p*_w - p*_e)A
   ```

3. **Solve pressure correction equation**
   ```
   Σ a_nb p'_nb - a_P p'_P = -Σ ṁ*_faces
   ```

4. **Correct pressure and velocity**
   ```
   p = p* + αp'
   u = u* + d_u(p'_w - p'_e)
   ```

5. **Update other properties**

6. **Check convergence**, return to step 2 if not converged

**Under-relaxation factors:**
- Pressure: α = 0.3-0.7
- Velocity: α = 0.7-0.9

### SIMPLEC (SIMPLE-Consistent)

**Modifications from SIMPLE:**
- Improved velocity correction formula
- Omits certain terms in denominator
- Allows higher under-relaxation factors
- Faster convergence

**Advantages:**
- More robust than SIMPLE
- Converges faster
- Can use α_p up to 0.9-1.0

### SIMPLER (SIMPLE-Revised)

**Key difference:**
- Solves pressure equation directly (not correction)
- Uses pressure correction only for velocity

**Steps:**
1. Solve momentum with guessed pressure
2. Solve pressure equation
3. Solve pressure correction for velocity
4. Update velocities

**Advantages:**
- More robust for some flows
- Better initial convergence
- Good when pressure field is important

### PISO (Pressure Implicit with Splitting of Operators)

Originally developed for transient flows, now used for steady-state.

**Steps:**
1. **Predictor:** Solve momentum with old pressure
2. **Corrector 1:** Solve pressure correction, update velocity
3. **Corrector 2:** Additional pressure correction
4. Can add more correctors if needed

**Advantages:**
- Better for transient problems
- Allows larger time steps
- Used in OpenFOAM
- Good for highly unsteady flows

**Disadvantages:**
- More expensive per iteration
- Multiple pressure corrections

### Fractional Step Method

Also called projection method or pressure correction method.

**Algorithm:**
1. **Momentum prediction** (without pressure gradient):
   ```
   u* = u^n + Δt[-(u·∇)u + ν∇²u + f]
   ```

2. **Pressure Poisson equation**:
   ```
   ∇²p^{n+1} = (1/Δt)∇·u*
   ```

3. **Velocity correction**:
   ```
   u^{n+1} = u* - Δt∇p^{n+1}
   ```

**Advantages:**
- Decouples velocity and pressure
- Efficient for incompressible flow
- Satisfies continuity exactly

## Compressible Flow Solution Methods

### Explicit Time Marching

**General form:**
```
U^{n+1} = U^n - (Δt/V)[R(U^n)]
```

Where R is the residual vector.

**Popular schemes:**

**Runge-Kutta Methods:**
- Multi-stage time integration
- 2nd, 3rd, or 4th order accuracy
- Good stability properties

**Example (4-stage RK):**
```
U^(1) = U^n - (α₁Δt/V)R^n
U^(2) = U^n - (α₂Δt/V)R^(1)
U^(3) = U^n - (α₃Δt/V)R^(2)
U^(n+1) = U^n - (α₄Δt/V)R^(3)
```

**Advantages:**
- Simple to implement
- Low memory requirements
- Easy to parallelize

**Disadvantages:**
- CFL restriction on time step
- Slow convergence to steady state
- Many time steps needed

### Implicit Time Marching

**General form:**
```
[I + (Δt/V)∂R/∂U]ΔU = -(Δt/V)R^n
```

**Solution approaches:**

**Approximate Factorization:**
- Factor the implicit operator
- Solve sequence of simpler systems
- AF1, AF2 schemes

**Newton-Raphson:**
- Full linearization
- Solve large sparse system
- Fast convergence but expensive

**Advantages:**
- Unconditionally stable
- Large time steps allowed
- Fast convergence to steady state

**Disadvantages:**
- Complex implementation
- Higher memory requirements
- Requires solving large systems

### Flux Splitting Methods

Decompose convective flux based on wave propagation.

**Flux Vector Splitting:**
- Split flux based on Jacobian eigenvalues
- Steger-Warming, van Leer schemes

**Flux Difference Splitting:**
- Based on Riemann problem
- Roe's scheme, HLL, HLLC

**Advantages:**
- Captures shock waves accurately
- Physically motivated
- Robust for compressible flows

## Multigrid Methods

Accelerate convergence by solving on multiple grid levels.

### V-Cycle Algorithm

```
procedure V-Cycle(level):
    if level == coarsest:
        Solve exactly
    else:
        Pre-smooth on current level
        Restrict residual to coarser level
        V-Cycle(level-1)
        Prolongate correction to finer level
        Post-smooth on current level
```

### FAS (Full Approximation Storage)

For non-linear problems:
- Transfer solution (not just correction) to coarse grid
- Solve non-linear equations on coarse grid
- Transfer correction back

**Advantages:**
- Dramatic speedup (10-100x)
- Grid-independent convergence
- Essential for large problems

**Implementation considerations:**
- Requires hierarchy of grids
- Grid transfer operators crucial
- Smoother selection important

## Convergence Acceleration Techniques

### Local Time Stepping

Use local CFL-limited time step for each cell:
```
Δt_i = CFL × (Δx_i/|u_i|)
```

**Advantages:**
- Faster convergence to steady state
- Adapts to local flow conditions

**Disadvantages:**
- Only for steady-state problems
- Time-accurate solution requires global time step

### Residual Smoothing

Smooth residual before updating solution:
```
R̃ = (1-ε)R + ε(R_neighbors)
```

**Advantages:**
- Allows larger time steps
- Improves stability
- Reduces high-frequency errors

### Enthalpy Damping

For high Mach number flows:
- Damps acoustic waves
- Allows larger time steps
- Improves robustness

### Preconditioning

Modify governing equations to improve conditioning.

**Low Mach number preconditioning:**
- Equalizes wave speeds
- Improves convergence at low Mach numbers
- Essential for all-speed flows

**Advantages:**
- Faster convergence
- Better for ill-conditioned problems
- Extends applicability

## Adaptive Methods

### Adaptive Time Stepping

Automatically adjust time step based on:
- Local CFL number
- Solution changes
- Convergence rate
- User-specified criteria

**Typical strategy:**
```
Δt^{n+1} = Δt^n × min(factor_max, max(factor_min, safety × (tol/error)^power))
```

### Adaptive Mesh Refinement (AMR)

Dynamically refine/coarsen mesh based on solution features.

**Refinement criteria:**
- Gradient-based (velocity, pressure, etc.)
- Error estimators
- Feature detection (shocks, vortices)

**Advantages:**
- Concentrates resolution where needed
- Reduces computational cost
- Improves accuracy

**Disadvantages:**
- Complex implementation
- Dynamic load balancing needed
- Mesh management overhead

## Parallel Computing

### Domain Decomposition

Split computational domain among processors.

**Approaches:**
- Geometric partitioning
- Graph-based partitioning (METIS)
- Load balancing

**Communication patterns:**
- Exchange boundary data
- MPI (Message Passing Interface)
- OpenMP for shared memory

**Scalability considerations:**
- Communication overhead
- Load balancing
- Partitioning quality

### GPU Acceleration

Utilize graphics processing units for CFD.

**Suitable operations:**
- Dense linear algebra
- Explicit time marching
- Matrix-vector products

**Challenges:**
- Memory bandwidth
- Limited memory
- Algorithm redesign needed

## Best Practices

### Solution Strategy Selection

**Incompressible flows:**
- Use pressure-based solver
- SIMPLE or PISO algorithm
- Segregated approach

**Compressible flows:**
- Use density-based solver
- Implicit time marching
- Coupled approach

**Transient flows:**
- PISO algorithm
- Appropriate time step
- Second-order time discretization

### Convergence Monitoring

Monitor:
- Residuals (all equations)
- Global imbalances (mass, energy)
- Solution variables at key locations
- Forces and moments (if relevant)

**Typical criteria:**
- Residuals < 10⁻³ to 10⁻⁶
- Imbalances < 0.1-1%
- Monitored variables stabilized

### Initialization

Good initial conditions accelerate convergence:
- Patch velocities based on physics
- Initialize pressure field
- Warm start from coarser mesh
- Use potential flow solution
