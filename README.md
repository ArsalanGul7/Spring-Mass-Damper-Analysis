# Spring-Mass-Damper Vibration Analysis (MATLAB)

This MATLAB project simulates and analyzes the behavior of a **1-DOF spring-mass-damper system** subjected to **harmonic excitation**. The script performs:

- Time-domain simulation using `ode45`
- Frequency-domain analysis using Fast Fourier Transform (FFT)
- Automatic interpretation of damping, resonance, and steady-state response
- Plots of displacement, velocity, and frequency spectrum
- Saves simulation output to `.mat` file

---

## System Overview

- **Mass (m):** 1.0 kg  
- **Damping Coefficient (c):** 2.0 Ns/m  
- **Stiffness (k):** 20.0 N/m  
- **Excitation Force:** `F(t) = 5·sin(3t)`  
- **Initial Conditions:** `[x(0) = 0,  ẋ(0) = 0]`  
- **Simulation Time:** 0 to 20 seconds

---

## Features

- Calculates natural frequency and damping ratio
- FFT-based frequency spectrum to find dominant frequency
- Interprets system as underdamped / overdamped / critically damped
- Detects resonance potential
- Saves all simulation variables into `SpringMassDamper_Analysis.mat`

---

## Plots

- **Displacement vs. Time**
- **Velocity vs. Time**
- **Frequency Spectrum of Displacement**

---

## How to Run

1. Clone or download this repository.
2. Open the file `Vibration_Analysis_of_a_Spring_Mass_Damper_System.m` in MATLAB.
3. Run the script.
4. View results in the MATLAB console and figures.
5. Output `.mat` file will be saved in the current folder.

---

## Files

| File                                                     | Description                                |
|----------------------------------------------------------|--------------------------------------------|
| `Vibration_Analysis_of_a_Spring_Mass_Damper_System.m`    | Main MATLAB script                         |
| `SpringMassDamper_Analysis.mat`                          | Output file with all variables             |
| `README.md`                                              | This readme file                           |
| `Displacement.jpg`                                       | Plot of Time-Domain Response: Displacement |
| `Velocity.jpg`                                           | Plot of Time-Domain Response: Velocity     |
| `Frequency Spectrum of Displacement.jpg`                 | Plot of Frequency Spectrum of Displacement |

---
