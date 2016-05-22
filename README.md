# 2D Discrete Inverse Spectral Problem

### Goal 
To calculate an approximate solution to the classic (Laplace-Beltrami) inverse spectral problem for discrete surfaces (genus 0 for now).

### Content
Included is a suite of MATLAB codes implementing the naive direct gradient descent approach.

test_script.m is the top-level script that generates results.

### People
Project envisioned, advised, and supervised by **Prof. Etienne Vouga** and **Prof. Keenan Crane**

Some codes here (on (sphere) meshes optimization and a demo of spherical harmonics) are not mine.

### Results
![ani#1](/i3_300_t2_abs(Y33(v))_e0.1-1p0.5.gif?raw=true "discrete Y33 spherical harmonic target with varying percent of eigenvalues used")
![ani#1](/i2_300_t2_abs(Y32(v))_e0.6p0.4-2.gif?raw=true "discrete Y32 spherical harmonic target with varying amount of deformation")

[//]: #![Result #1](/rand_0.512.png?raw=true "target constructed with small random discrete conformal factors")
[//]: #![Result #2](/Y10_0.512.png?raw=true "discrete Y10 spherical harmonic target, wat")
[//]: #![Result #3](/Y10_0.125.png?raw=true "discrete Y10 spherical harmonic target (with smaller perturbation, a.k.a. ellipsoid)")
[//]: #![Result #4](/Y20_0.512.png?raw=true "discrete Y20 spherical harmonic target, converged at a local minimum")
[//]: #![Result #5](/Y20_0.216.png?raw=true "discrete Y20 spherical harmonic target (with less perturbation)")
[//]: #![Result #6](/Y33_0.343.png?raw=true "discrete Y33 spherical harmonic target, 500 vtx mesh result")