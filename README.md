# Kinetic Exclusion Process & Macroscopic Fluctuation Theory (MFT)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![Fortran](https://img.shields.io/badge/Fortran-gfortran-purple.svg)](https://gcc.gnu.org/fortran/)

This repository contains the source code, numerical solvers, and raw data associated with my Bachelor's Thesis (Trabajo Fin de Grado) in the Double Degree in Physics and Mathematics at the University of Granada (UGR).

**Author:** Sharan Haresh Belani  
**Tutor:** Dr. Pablo Ignacio Hurtado Fernández  

## Overview

This project studies the non-equilibrium transport properties and macroscopic fluctuations of the **Kinetic Exclusion Process (KEP)**. By dynamically coupling density and energy fields through microscopic jump rules, the KEP captures complex cross-transport phenomena such as the Soret and Dufour effects. 

The repository includes Kinetic Monte Carlo (KMC) simulations written in Fortran, macroscopic hydrodynamic solvers, and Python visualization scripts to compare simulation data with theoretical predictions.

## Repository Structure

The repository is organized into folders based on the parameter sweeps and physical studies performed during the research:

*  **`Caso_beta_0/`, `Caso_beta_05/`, `Caso_beta_1/`**
    * Systematic studies of the collision parameter $\beta$ (0.0, 0.5, and 1.0). 
    * Each folder contains subdirectories for different temperature gradients and density boundary conditions.
    * Includes the KMC simulation code (`kep_simple_fixed_...f90`), hydrodynamic solvers (`hydro_solver...f90`), Python plotting scripts, and the resulting density/temperature data (`.dat` and `.csv`).
*  **`Estudio_densidad/`**
    * Parametric study isolating the effect of the average background density ($\rho = 0.2, 0.5, 0.9$).
    * Contains comparison plots (`comparacion_density.png`) and the corresponding source codes.
*  **`Estudio_termico/`**
    * Parametric study isolating the effect of extreme temperature gradients ($T = 2, 5, 10$).
    * Contains comparison plots (`comparacion_thermal.png`).
*  **`KMP_Validation/` & `SSEP_Validation/`**
    * Validation of the KMC algorithm against exactly solvable canonical models.
    * Contains `kmp_sim.f90` (pure heat conduction limit) and `ssep_sim.f90` (pure mass diffusion limit), along with validation plots.
*  **`Caso_Canonico/`**
    * Study of the canonical ensemble where the Local Equilibrium Approximation (LEA) breaks down, including the breakdown analysis plots (`analisis_ruptura_LEA.png`).
*  **`FiniteSizeEffects/` & `Variable_beta/`**
    * Additional data and scripts analyzing boundary interactions and continuous variations of the coupling parameter.
