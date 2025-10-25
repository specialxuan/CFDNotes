# Discretization Methods

## Overview

Discretization is the process of converting continuous partial differential equations (PDEs) into discrete algebraic equations that can be solved numerically on a computer. The domain is divided into small control volumes or elements, and the governing equations are approximated at discrete locations.

## Finite Difference Method (FDM)

The Finite Difference Method approximates derivatives using Taylor series expansions on structured grids.

### Forward Difference
```
∂f/∂x ≈ (f_{i+1} - f_i)/Δx
```

### Backward Difference
```
∂f/∂x ≈ (f_i - f_{i-1})/Δx
```

### Central Difference
```
∂f/∂x ≈ (f_{i+1} - f_{i-1})/(2Δx)
```

### Second Derivative
```
∂²f/∂x² ≈ (f_{i+1} - 2f_i + f_{i-1})/Δx²
```

### Advantages
- Simple to implement
- High accuracy on structured grids
- Straightforward error analysis

### Disadvantages
- Limited to simple geometries
- Difficult to apply on irregular domains
- Requires structured meshes

## Finite Volume Method (FVM)

The Finite Volume Method divides the domain into control volumes and applies conservation laws in integral form.

### Discretization Process

1. Divide domain into control volumes
2. Integrate governing equations over each control volume
3. Apply divergence theorem to convert volume integrals to surface integrals
4. Approximate fluxes at cell faces
5. Assemble algebraic equations

### General Form

For a generic transport equation:
```
∂(ρφ)/∂t + ∇·(ρuφ) = ∇·(Γ∇φ) + S
```

Integrated over control volume:
```
∫_V ∂(ρφ)/∂t dV + ∫_V ∇·(ρuφ) dV = ∫_V ∇·(Γ∇φ) dV + ∫_V S dV
```

Using divergence theorem:
```
∂(ρφV)/∂t + Σ_faces (ρuφA)_f = Σ_faces (Γ∇φ·A)_f + SV
```

### Advantages
- Naturally conservative (conserves mass, momentum, energy)
- Can handle complex geometries
- Works on structured and unstructured grids
- Most widely used in commercial CFD codes

### Disadvantages
- Lower accuracy compared to FEM on irregular grids
- Complex implementation for higher-order schemes

## Finite Element Method (FEM)

The Finite Element Method approximates the solution using piecewise polynomial functions over elements.

### Discretization Process

1. Divide domain into elements
2. Choose basis functions (shape functions)
3. Apply weighted residual method (Galerkin method)
4. Assemble element matrices into global system
5. Apply boundary conditions
6. Solve resulting system of equations

### Weak Form

Starting from strong form:
```
L(φ) = f
```

Multiply by test function w and integrate:
```
∫_Ω w·L(φ) dΩ = ∫_Ω w·f dΩ
```

Apply integration by parts to reduce derivative order.

### Shape Functions

Common element types:
- **Linear elements**: First-order polynomials
- **Quadratic elements**: Second-order polynomials
- **Triangular elements**: 2D domains
- **Tetrahedral elements**: 3D domains
- **Quadrilateral elements**: 2D structured
- **Hexahedral elements**: 3D structured

### Advantages
- High accuracy with adaptive refinement
- Flexible for complex geometries
- Strong mathematical foundation
- Excellent for structural mechanics

### Disadvantages
- More complex to implement
- Higher computational cost
- Not inherently conservative (requires special treatment)
- Less commonly used for fluid dynamics

## Spectral Methods

Spectral methods use global basis functions (e.g., Fourier series, Chebyshev polynomials) to represent the solution.

### Fourier Spectral Method

Approximate solution as:
```
φ(x) = Σ_k â_k e^{ikx}
```

### Chebyshev Spectral Method

Approximate solution using Chebyshev polynomials:
```
φ(x) = Σ_n a_n T_n(x)
```

### Advantages
- Exponential convergence for smooth solutions
- Very high accuracy
- Efficient for periodic problems

### Disadvantages
- Limited to simple geometries
- Poor performance for discontinuous solutions
- Requires smooth solutions for optimal performance

## Mesh Types

### Structured Meshes
- Grid lines form a regular pattern
- Can be mapped to Cartesian grid
- Efficient storage and computation
- Limited to simple geometries

### Unstructured Meshes
- Arbitrary connectivity between cells
- Flexible for complex geometries
- Higher memory requirements
- More complex data structures

### Hybrid Meshes
- Combination of structured and unstructured regions
- Structured mesh near walls (boundary layers)
- Unstructured mesh in complex regions

## Time Discretization

### Explicit Methods
```
φ^{n+1} = φ^n + Δt·F(φ^n)
```

- Simple to implement
- Conditionally stable (CFL condition)
- Small time steps required

### Implicit Methods
```
φ^{n+1} = φ^n + Δt·F(φ^{n+1})
```

- Unconditionally stable
- Larger time steps allowed
- Requires solving system of equations

### Crank-Nicolson (Semi-Implicit)
```
φ^{n+1} = φ^n + Δt·[θF(φ^{n+1}) + (1-θ)F(φ^n)]
```

- Second-order accurate when θ = 0.5
- Good stability properties
- Widely used for transient simulations

## Accuracy and Convergence

### Truncation Error

The difference between the exact derivative and its discrete approximation:
```
Truncation Error = O(Δx^p)
```

Where p is the order of accuracy.

### Grid Convergence

As mesh is refined (Δx → 0):
- Solution should converge to exact solution
- Monitor solution changes with refinement
- Use Richardson extrapolation

### CFL Condition

For explicit time marching:
```
CFL = |u|Δt/Δx ≤ CFL_max
```

Typical CFL_max ≈ 1 for stability.
