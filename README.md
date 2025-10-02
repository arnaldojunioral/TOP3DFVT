# ✨ Top3DFVT

**Topology Optimization of 3D Elastic Structures employing the Finite-Volume Theory (FVT)**

This repository provides a **MATLAB implementation** of topology optimization for 3D continuum elastic structures employing the **Finite-Volume Theory (FVT)**. The algorithm optimizes the material distribution in a discretized domain to achieve a minimum structural compliance.

Supported material interpolation methods:
- **SIMP** (Solid Isotropic Material with Penalization)  
- **RAMP** (Rational Approximation of Material Properties)

or
 
- **GSS** (Gray Scale Suppression - Voigt model) - strategy to penalize intermediate densities, effectively guiding the optimization toward crisp, manufacturable topologies.

<p align="center">
  <img width="285" height="268" alt="image" src="https://github.com/user-attachments/assets/f57dceaa-2773-4d08-992f-6858ca71c8fa" /><br>
  <em>Metaphorical illustration of human evolution and topology optimization with Gray Scale Suppression (GSS)</em>
</p>

----

## 📌 Features

- 3D topology optimization of elastic structures.
- It supports either a continuation scheme applied to penalization factors or a fixed penalization approach.
- No filtering solutions.

- Multiple **benchmark problems**:
  - Cantilever beam
  - MBB beam
  - Bridge

- 3D optimized design.
- Black-and-white fraction calculation and XY-plane symmetry check.

----

- ## ⚙️ Requirements

The implementation is fully compatible with both **MATLAB (R2015 or later)** and **GNU Octave (version 6.4 or later)**. No additional toolboxes are required, ensuring that the code can be executed in a standard installation of either environment.

----

## 🚀 Getting started

Save the [TOP3DFVT.m](https://raw.githubusercontent.com/arnaldojunioral/TOP3DFVT/main/TOP3DFVT.m) program (16.7 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

**Top3DFVT(n1, n2, n3, volfrac, model)**

where **n1**, **n2**, and **n3** define the number of subvolumes along the x<sub>1</sub>, x<sub>2</sub>, and x<sub>3</sub> directions, respectively. These parameters specify the discretization of the three-dimensional domain, as illustrated in the figure below. **volfrac** is the volume fraction constraint of solid material, and **model** defines the penalization strategy: **SIMP**, **RAMP**, or **GSS**.

<!-- <p align="center">
<img width="350" height="350" alt="image" src="https://github.com/user-attachments/assets/3d92838e-2fcb-40f7-b0da-80d891ec62d6" />
</p> -->
<p align="center">
<img width="510" height="280" alt="image" src="https://github.com/user-attachments/assets/9efeaf36-5ae0-45d4-b3fc-7a83c7a8b952" />
</p>

The table below summarizes the key input parameters used in the simulation, including beam geometry, material properties, loading conditions, and visualization settings.

#### Model Parameters

| Parameter       | Description                                          | Unit             |
|-----------------|------------------------------------------------------|------------------|
| L             | Beam length                                          | mm               |
| H             | Beam height                                          | mm               |
| B             | Beam width                                           | mm               |
| E             | Young's modulus (material stiffness)                 | MPa              |
| nu            | Poisson's ratio                                      | –                |
| P             | Applied load (negative indicates downward force)     | N                |
| eta          | 1/3      | Damping factor                                    |
| tol          | 1e-3     | Convergence tolerance                                       |
| maxit        | 1000     | Maximum number of iterations                                |
| threshold    | 0.01     | Density threshold for display                               |
| displayflag  | true     | Display during optimization (true = live updates)         |
| pb           | Cantilever | Problem type                                              |

----

## 📌 Example

In the following example, the model consists of a structured mesh discretized into 30 × 10 × 3  subvolumes, with a solid material volume fraction constrained to 50%.
The GSS strategy was employed.

<p align="center">
<img width="500" height="280" alt="image" src="https://github.com/user-attachments/assets/3286338d-7d7f-4e84-8526-e4e316ec63df" />
</p>

Run the main function:

```matlab

Top3DFVT(30, 10, 0.5, 'GSS');  

[L, H, B] = deal(120, 40, 20);   % Beam dimensions
pb = 'Cantilever'; 
E0 = 1;        % Young's modulus of solid material
Emin = 1e-9;   % Young's modulus of void (numerical stability)
nu = 0.3;      % Poisson ratio
P = -1;        % Applied load
eta = 1/3;     % Filtering damping factor
tol = 1e-3;    % Convergence tolerance
maxit = 1000;  % Maximum iterations
threshold = 0.01; % Density threshold for display
displayflag = true; % Live plot flag

````

----

## 🎥 Topology evolution

<p align="center">
<b> Gray Scale Suppression<b> </b>
</p>
 
<p align="center">
  <img width="500" height="300" alt="image" src="https://github.com/user-attachments/assets/16f8fb84-74f8-44b9-b356-09a90b7796d1" /><br>
  B&W: 0.999<br>
  Symmetry across x<sub>1</sub>x<sub>2</sub>-plane: true
</p>

----
