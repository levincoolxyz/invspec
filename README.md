# 2D Discrete Inverse Spectral Problem

### Goal 
To calculate an approximate solution to the classic (Laplace-Beltrami) inverse spectral problem for discrete surfaces (genus 0 for now).

### Content
Included is a suite of MATLAB codes implementing the naive direct gradient descent approach.
test_script.m is the top-level script that generates results.

### People
Project envisioned, advised, and supervised by **Prof. Etienne Vouga** and **Prof. Keenan Crane**

Some codes here (on (sphere) meshes optimization and a demo of spherical harmonics) are not mine.

### Some Results
![Result #1](/rand_0.8.png?raw=true "target constructed with small random discrete conformal factors")
![Result #3](/Y10_0.8.png?raw=true "discrete Y10 spherical harmonic target, wat")
![Result #3](/Y10_0.5.png?raw=true "discrete Y10 spherical harmonic target (with smaller perturbation, a.k.a. ellipsoid)")
![Result #2](/Y20_0.8.png?raw=true "discrete Y20 spherical harmonic target, converged at a local minimum")
![Result #3](/Y20_0.6.png?raw=true "discrete Y20 spherical harmonic target (with less perturbation)")
![Result #3](/Y33_0.8.png?raw=true "discrete Y33 spherical harmonic target")
![Result #3](/Y33_0.7.png?raw=true "discrete Y33 spherical harmonic target, finer mesh less perturbation")