# 2D Truss Analysis using MATLAB

This repository contains a complete MATLAB script for the static analysis of 2D truss structures using the **Direct Stiffness Method**.

## Features

- Input any number of nodes and members.
- Automatic assembly of the global stiffness matrix.
- Support for fixed and roller boundary conditions.
- Calculation of nodal displacements (Ux, Uy).
- Calculation of reaction forces at supports.
- Calculation of internal axial forces and normal stresses for each member.
- Stability check (detects singular stiffness matrix).

## How to Use

1. Clone or download the `Truss_Analysis.m` file.
2. Open the file in MATLAB (version R2020a or later recommended).
3. Modify the **"Input Data"** section at the top of the script to match your structure:
   - `X`, `Y`: Coordinates of nodes.
   - `Conn`: Connectivity matrix (member connections).
   - `E`, `A`: Young's modulus and cross-sectional area for each member.
   - `F_ext`: External forces on nodes.
   - `BC`: Boundary conditions (`1` = fixed, `0` = free).
4. Run the script.
5. View the results in the MATLAB Command Window.

## Input Format Example

```matlab
X = [0; 4; 6; 2; 3];      % Node X coordinates
Y = [0; 0; 0; 3; 5];      % Node Y coordinates

Conn = [1 2; 2 3; 1 4; 2 4; 3 4; 4 5; 3 5];  % Connections

E = [200e9; 200e9; 70e9; ...];   % Young's modulus (Pa)
A = [0.001; 0.001; 0.002; ...];  % Area (m^2)

F_ext = [0 0; 0 0; ...];         % External forces [Fx, Fy] (N)
BC = [1 1; 1 1; 0 0; ...];       % Boundary conditions
