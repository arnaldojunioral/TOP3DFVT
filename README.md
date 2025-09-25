# Top3DFVT

**Topology Optimization of 3D Elastic Structures employing the Finite-Volume Theory (FVT)**

This repository provides a **MATLAB/Octave implementation** of topology optimization for 3D continuum elastic structures employing the **Finite-Volume Theory (FVT)**. The algorithm optimizes the material distribution in a discretized domain to achieve desired structural properties.

---

## ⚙️ Features

- 3D topology optimization of elastic structures.
- Supports different **interpolation models**:
  - **SIMP** (Solid Isotropic Material with Penalization)
  - **RAMP** (Rational Approximation of Material Properties)
  - **GSS** (Grayscale Suppression - linear interpolation)

- Multiple **benchmark problems**:
  - Cantilever beam
  - MBB beam
  - Bridge

- Black-and-white fraction calculation and XY-plane symmetry check.
- Live 3D visualization during optimization.
- Output: optimized density matrix (`Density.mat`) and optimized design.

---

## 📌 User-Defined Parameters

```matlab
nx = 120;      % number of subvolumes in x-direction
ny = 40;       % number of subvolumes in y-direction
nz = 20;       % number of subvolumes in z-direction
volfrac = 0.5; % target solid material volume fraction

model = 'GSS'; % 'SIMP', 'RAMP', or 'GSS'
pb = 'Cantilever'; % Problem type: 'Cantilever', 'MBB', 'Bridge'

E0 = 1;        % Young's modulus of solid material
Emin = 1e-9;   % Young's modulus of void (numerical stability)
nu = 0.3;      % Poisson ratio
P = -1;        % Applied load
eta = 1/3;     % Filtering damping factor
tol = 1e-3;    % Convergence tolerance
maxit = 1000;  % Maximum iterations
threshold = 0.01; % Density threshold for display
displayflag = true; % Live plot flag
