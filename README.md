# Code for Paper "Advancing 3D Line Registration with Geometric Algebra: A Framework for Unified Transformations"

## Overview
This MATLAB project implements a unified framework for 3D line registration by using Geometric Algebra (GA) techniques combined with standard SVD methods. The script maps 3D Euclidean line data into a 4D spherical space using the 1DUp approach within Conformal Geometric Algerba (CGA), applies transformations using GA and recovers the 3D rotation and translation components. It then compares the GA-based registration with a conventional SVD point-based approach.

## Features
- **3D Line Generation:** Generates multiple 3D lines with configurable numbers of lines and points per line.
- **Noise Addition:** Optionally introduces Gaussian noise to simulate real-world measurement errors.
- **TLS Fitting:** Uses Total Least Squares (TLS) to fit lines to noisy data.
- **Spherical Mapping:** Transforms 3D Euclidean points into 4D spherical space.
- **Geometric Algebra Processing:** Applies GA-based transformations and recovers the motor (transformation) that encapsulates both rotation and translation.
- **SVD Comparison:** Compares the GA-based registration with a standard SVD-based rigid transformation.
- **Error Analysis:** Computes rotational and translational errors for performance evaluation.
- **Modular Design:** Organized in a functional style with pure helper functions, making the code easy to read, test and maintain.

## Requirements
- MATLAB (R2020a or later is recommended)
- Clifford Algebra Toolbox (ensure it is installed and added to your MATLAB path)

## How to Run
1. **Install Requirements:**  
   Make sure you have MATLAB and the Clifford Algebra Toolbox installed.

2. **Save the Script:**  
   Save the provided code into a single file named `line_registration.m`.

3. **Run the Script:**  
   Open MATLAB, navigate to the directory containing `line_registration.m`, and execute:
   ```matlab
   line_registration
   ```
   The script will:
   - Generate multiple configurations of 3D lines.
   - Optionally add noise and perform TLS fitting.
   - Map the points to 4D spherical space and apply GA-based transformations.
   - Recover the transformation parameters (rotation and translation).
   - Compare these with the results from a standard SVD-based approach.
   - Plot the error metrics for both methods.

## File Structure
All the code is contained within a single file (`line_registration.m`). The script is organized as follows:
- **Initial Setup & Parameter Definition:**  
  Sets up simulation parameters, initializes the Clifford Algebra, and preallocates result structures.
- **Main Loop:**  
  Iterates over different configurations (varying numbers of lines and points per line) to generate lines, add noise, perform TLS fitting, map to spherical space, recover transformations and compute errors.
- **Error Calculation & Plotting:**  
  Calculates and plots rotational and translational error metrics for both the GA-based and SVD-based methods.
- **Helper Functions:**  
  Contains all supporting functions for line generation, noise addition, spherical mapping, transformation recovery and error computation. These functions follow a pure-functional style, meaning they take inputs and return outputs without side effects.

## Contact

For any questions or issues with the code, please contact:  
**[Haris Matsantonis]**  
**[cm2128@cam.ac.uk]**


