# Turbulence Modeling

## Overview

Turbulent flows are characterized by chaotic, random, three-dimensional vortical motion occurring across a wide range of spatial and temporal scales. Most engineering flows are turbulent, making turbulence modeling essential for practical CFD applications.

## Characteristics of Turbulent Flow

### Key Features
- **Irregular and chaotic**: Random fluctuations in velocity, pressure, and other flow properties
- **Diffusive**: Enhanced mixing and transport compared to laminar flow
- **Three-dimensional**: Vortical structures in all three dimensions
- **Wide range of scales**: From large eddies to small dissipative scales
- **High Reynolds number**: Typically Re > 2300 for pipe flow

### Turbulent Kinetic Energy
```
k = 1/2(u'² + v'² + w'²)
```
Where u', v', w' are velocity fluctuations.

### Energy Cascade
- Large eddies extract energy from mean flow
- Energy cascades to smaller eddies
- Smallest eddies dissipate energy as heat
- Kolmogorov scale: smallest dissipative eddies

## Approaches to Turbulence Simulation

### Direct Numerical Simulation (DNS)

Solves the Navier-Stokes equations without any turbulence modeling.

#### Characteristics
- Resolves all turbulent scales from largest to smallest
- No modeling assumptions
- Provides complete flow information
- Gold standard for turbulence research

#### Requirements
Grid resolution must satisfy:
```
Δx ≈ η = (ν³/ε)^(1/4)  (Kolmogorov scale)
N ≈ Re^(9/4)  (Number of grid points)
```

#### Limitations
- Extremely expensive: Re⁹/⁴ computational cost
- Limited to low Reynolds numbers and simple geometries
- Impractical for engineering applications
- Requires massive computational resources

### Large Eddy Simulation (LES)

Resolves large-scale turbulent motions and models small-scale (subgrid-scale) turbulence.

#### Filtering
```
φ̄(x) = ∫ φ(x')G(x-x', Δ)dx'
```
Where G is a filter function and Δ is the filter width.

#### Subgrid-Scale (SGS) Models

**Smagorinsky Model**
```
τ_ij - 1/3τ_kk δ_ij = -2ν_t S̄_ij
ν_t = (C_s Δ)² |S̄|
```
Where C_s ≈ 0.1-0.2.

**Dynamic Smagorinsky Model**
- Computes model coefficient dynamically
- More accurate than standard Smagorinsky
- Adjusts to local flow conditions

**WALE Model (Wall-Adapting Local Eddy)**
- Better near-wall behavior
- Accounts for strain and rotation
- No explicit damping functions needed

#### Advantages
- Resolves energy-containing scales
- More accurate than RANS
- Suitable for unsteady, separated flows
- Captures large-scale turbulent structures

#### Disadvantages
- Still computationally expensive
- Requires fine mesh near walls
- Sensitive to boundary conditions
- Time-dependent (long averaging times needed)

### Reynolds-Averaged Navier-Stokes (RANS)

Solves time-averaged equations with turbulence models for all scales.

#### Reynolds Decomposition
```
φ = φ̄ + φ'
```
Where φ̄ is the mean and φ' is the fluctuation.

#### Reynolds Stresses

After averaging, additional terms appear:
```
τ_ij = -ρu'_i u'_j  (Reynolds stress tensor)
```

This creates the closure problem: more unknowns than equations.

## RANS Turbulence Models

### Zero-Equation Models

#### Mixing Length Model
```
ν_t = l_m² |∂u/∂y|
```

- Simplest turbulence model
- Algebraic (no transport equations)
- Limited applicability
- Requires prescription of mixing length l_m

### One-Equation Models

#### Spalart-Allmaras Model

Solves transport equation for modified eddy viscosity ν̃:
```
∂ν̃/∂t + u_j ∂ν̃/∂x_j = Production - Destruction + Diffusion
```

**Characteristics:**
- Single equation
- Economical
- Good for aerospace applications
- Works well for boundary layers
- Limited for separated flows

**Advantages:**
- Robust and stable
- Low computational cost
- Good for attached flows

**Disadvantages:**
- Limited for complex flows
- Not ideal for free shear flows
- Single-scale model

### Two-Equation Models

Most widely used turbulence models in engineering.

#### k-ε Model (Standard)

Transport equations for turbulent kinetic energy k and dissipation rate ε:

**k-equation:**
```
∂k/∂t + u_j ∂k/∂x_j = P_k - ε + ∂/∂x_j[(ν + ν_t/σ_k) ∂k/∂x_j]
```

**ε-equation:**
```
∂ε/∂t + u_j ∂ε/∂x_j = C_1ε (ε/k)P_k - C_2ε (ε²/k) + ∂/∂x_j[(ν + ν_t/σ_ε) ∂ε/∂x_j]
```

**Eddy viscosity:**
```
ν_t = C_μ k²/ε
```

**Model constants:**
```
C_μ = 0.09, C_1ε = 1.44, C_2ε = 1.92, σ_k = 1.0, σ_ε = 1.3
```

**Advantages:**
- Well-established and robust
- Good for fully turbulent flows
- Widely validated

**Disadvantages:**
- Poor near walls (requires wall functions)
- Overpredicts turbulence in stagnation regions
- Limited for swirling and rotating flows

#### k-ε Model Variants

**RNG k-ε:**
- Derived using Renormalization Group theory
- Better for swirling flows
- Improved accuracy for rapidly strained flows

**Realizable k-ε:**
- Ensures realizability of Reynolds stresses
- Better for separated flows and jets
- Modified dissipation equation

#### k-ω Model (Standard)

Uses specific dissipation rate ω = ε/k instead of ε:

**k-equation:**
```
∂k/∂t + u_j ∂k/∂x_j = P_k - β*kω + ∂/∂x_j[(ν + ν_t/σ_k) ∂k/∂x_j]
```

**ω-equation:**
```
∂ω/∂t + u_j ∂ω/∂x_j = α(ω/k)P_k - βω² + ∂/∂x_j[(ν + ν_t/σ_ω) ∂ω/∂x_j]
```

**Advantages:**
- Accurate near walls
- Good for adverse pressure gradients
- No wall functions needed

**Disadvantages:**
- Sensitive to freestream values of ω
- Less accurate in free shear layers

#### SST k-ω Model (Shear Stress Transport)

Combines k-ω near walls and k-ε in free stream:

**Key Features:**
- Blending function switches between models
- Limits eddy viscosity in adverse pressure gradients
- Industry standard for many applications

**Advantages:**
- Best of both k-ω and k-ε
- Excellent near-wall behavior
- Good for separated flows
- Widely used and validated

**Disadvantages:**
- More complex than standard models
- Slightly higher computational cost

### Reynolds Stress Models (RSM)

Solve transport equations for each Reynolds stress component.

#### Seven Transport Equations
- Six Reynolds stress components: u'_i u'_j
- One equation for dissipation rate ε or ω

**Advantages:**
- Most physically accurate RANS approach
- Accounts for stress anisotropy
- Better for complex flows (swirl, rotation, curvature)
- No assumption of eddy viscosity

**Disadvantages:**
- Computationally expensive (7 equations)
- Less robust than two-equation models
- Requires careful numerical treatment
- More complex to implement

## Wall Treatment

### Wall Functions

Bridge the near-wall region without resolving it:

**Law of the Wall:**
```
u⁺ = y⁺  (viscous sublayer, y⁺ < 5)
u⁺ = 1/κ ln(y⁺) + B  (log layer, 30 < y⁺ < 300)
```

Where:
```
y⁺ = yuτ/ν
u⁺ = u/uτ
uτ = √(τw/ρ)  (friction velocity)
```

**Standard Wall Functions:**
- First grid point at 30 < y⁺ < 300
- Economical (coarse near-wall mesh)
- Good for attached flows
- Less accurate for complex flows

**Enhanced Wall Functions:**
- Blend between wall functions and low-Re formulation
- More robust
- Work for wider range of y⁺

### Low-Reynolds Number Models

Resolve flow down to the wall:

**Requirements:**
- First grid point at y⁺ ≈ 1
- Multiple points in viscous sublayer (y⁺ < 5)
- Fine mesh near walls

**Advantages:**
- More accurate near walls
- Better for separated flows
- No wall function assumptions

**Disadvantages:**
- Requires very fine mesh
- Higher computational cost
- Grid sensitivity

## Transition Modeling

Predicting laminar-to-turbulent transition:

### Transition Models

**γ-Reθ Model:**
- Two additional transport equations
- Predicts transition location
- Works with SST k-ω

**Applications:**
- Low Reynolds number airfoils
- Turbomachinery
- Heat transfer prediction

## Best Practices

### Model Selection Guidelines

- **General purpose:** SST k-ω
- **External aerodynamics:** Spalart-Allmaras or SST k-ω
- **Internal flows:** k-ε or SST k-ω
- **Complex flows (swirl, rotation):** RSM
- **High accuracy needed:** LES
- **Transition important:** γ-Reθ transition model

### Grid Requirements

- k-ε with wall functions: y⁺ = 30-300
- k-ω or SST: y⁺ ≈ 1
- LES: Δx⁺ < 100, Δz⁺ < 20, y⁺ < 1

### Validation

- Compare with experimental data
- Check grid independence
- Verify realizability constraints
- Monitor turbulence quantities
