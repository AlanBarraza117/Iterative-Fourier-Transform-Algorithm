# Iterative Fourier Transform Algorithm (IFTA)

MATLAB implementation of the Iterative Fourier Transform Algorithm for designing Diffractive Optical Elements (DOE).

## Description

The algorithm reconstructs a target intensity pattern using iterative Fourier transforms. It alternates between the DOE plane and the output plane, constraining the amplitude at each step while updating the phase.

## Algorithm Steps

1. Apply FFT from the DOE plane to the output plane
2. Substitute the amplitude of |U_out| with the target amplitude, preserving the phase
3. Apply IFFT to go back to the DOE plane
4. Substitute the amplitude from |U_back| with the original amplitude (preserve phase)
5. Repeat until convergence

## Metrics

- **eta** - Diffraction Efficiency: Fraction of power reaching the signal window
- **UE** - Uniformity Error: How much the results differ from the target

## Usage

Run `IFTA.m` in MATLAB.

## Author

Jose Alan Barraza Villaverde - iPSRS (Machine Learning and Photonics)
