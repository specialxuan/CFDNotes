# Governing Equations

## Overview

The governing equations of fluid dynamics are derived from the fundamental physical principles of conservation of mass, momentum, and energy. These equations form the basis of all CFD simulations.

## Conservation of Mass (Continuity Equation)

The continuity equation expresses the conservation of mass in a fluid flow.

### General Form

```
∂ρ/∂t + ∇·(ρu) = 0
```

Where:
- ρ = density
- t = time
- u = velocity vector
- ∇· = divergence operator

### Incompressible Flow

For incompressible flows (constant density):

```
∇·u = 0
```

This simplifies to:
```
∂u/∂x + ∂v/∂y + ∂w/∂z = 0
```

## Conservation of Momentum (Navier-Stokes Equations)

The Navier-Stokes equations describe the motion of viscous fluid substances.

### General Form

```
ρ(∂u/∂t + u·∇u) = -∇p + ∇·τ + f
```

Where:
- p = pressure
- τ = viscous stress tensor
- f = body forces (e.g., gravity)

### Expanded Form (Incompressible, Newtonian Fluid)

```
ρ(∂u/∂t + u·∇u) = -∇p + μ∇²u + f
```

Where:
- μ = dynamic viscosity
- ∇² = Laplacian operator

### Component Form (3D Cartesian)

**X-momentum:**
```
ρ(∂u/∂t + u∂u/∂x + v∂u/∂y + w∂u/∂z) = -∂p/∂x + μ(∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²) + f_x
```

**Y-momentum:**
```
ρ(∂v/∂t + u∂v/∂x + v∂v/∂y + w∂v/∂z) = -∂p/∂y + μ(∂²v/∂x² + ∂²v/∂y² + ∂²v/∂z²) + f_y
```

**Z-momentum:**
```
ρ(∂w/∂t + u∂w/∂x + v∂w/∂y + w∂w/∂z) = -∂p/∂z + μ(∂²w/∂x² + ∂²w/∂y² + ∂²w/∂z²) + f_z
```

## Conservation of Energy

The energy equation accounts for heat transfer and work done by the fluid.

### General Form

```
ρc_p(∂T/∂t + u·∇T) = ∇·(k∇T) + Φ + Q
```

Where:
- T = temperature
- c_p = specific heat at constant pressure
- k = thermal conductivity
- Φ = viscous dissipation
- Q = heat source/sink

## Simplified Equations

### Euler Equations

For inviscid (frictionless) flows, the Navier-Stokes equations reduce to the Euler equations:

```
∂ρ/∂t + ∇·(ρu) = 0
ρ(∂u/∂t + u·∇u) = -∇p + f
```

### Stokes Flow

For very slow flows (creeping flow), inertial terms can be neglected:

```
0 = -∇p + μ∇²u + f
```

### Potential Flow

For irrotational, inviscid flows:

```
∇²φ = 0
```

Where φ is the velocity potential (u = ∇φ).

## Dimensionless Numbers

Important dimensionless parameters that characterize fluid flows:

### Reynolds Number (Re)
```
Re = ρUL/μ = UL/ν
```
- Ratio of inertial to viscous forces
- Indicates laminar vs. turbulent flow

### Mach Number (Ma)
```
Ma = U/c
```
- Ratio of flow velocity to speed of sound
- Classifies compressibility effects

### Froude Number (Fr)
```
Fr = U/√(gL)
```
- Ratio of inertial to gravitational forces
- Important for free surface flows

### Prandtl Number (Pr)
```
Pr = ν/α = μc_p/k
```
- Ratio of momentum to thermal diffusivity
- Characterizes heat transfer

### Grashof Number (Gr)
```
Gr = gβΔTL³/ν²
```
- Ratio of buoyancy to viscous forces
- Important for natural convection

## Boundary Conditions

Common boundary conditions used in CFD:

1. **Velocity Inlet**: Specify velocity components
2. **Pressure Inlet/Outlet**: Specify pressure
3. **Wall**: No-slip (u=0) or slip conditions
4. **Symmetry**: Zero normal velocity and gradients
5. **Periodic**: Flow repeats spatially
6. **Free Surface**: Pressure equals atmospheric

## Initial Conditions

For transient simulations, initial conditions must be specified:
- Velocity field at t=0
- Pressure field at t=0
- Temperature field at t=0 (if energy equation is solved)
