# 2D Discrete Inverse Spectral Problem

### Goal 
To calculate an approximate solution to the classic (Laplace-Beltrami) inverse spectral problem for discrete surfaces (genus 0 for now).

### Content
Included is a suite of MATLAB codes implementing the naive direct gradient descent approach.

test_script.m is the top-level script that generates results.

### People
Project envisioned, advised, and supervised by **Prof. Etienne Vouga** and **Prof. Keenan Crane**

Some codes here (on (sphere) meshes optimization and a demo of spherical harmonics) are not mine.

### Procedures (tentative)
1. (conformalized) mean curvature flow of target mesh onto "spherical" mesh with a target set of conformal factors

to do: rescale volume? spectral convergence?

2. descent search for some conformal factors that optimize spectral similarity

3. fit conformal factors on sphere back to a resulting mesh that optimize edge lengths targets from the factors

### Results
Aiming for cMCF spectrum

![ani#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/i2_300_t2_abs(Y33(v))_e0.1-1p0.5.gif "discrete Y33 spherical harmonic target with varying percent of eigenvalues used")
![ani#2](https://raw.githubusercontent.com/levincoolxyz/invspec/master/i2_300_t2_abs(Y32(v))_e0.1p0.5-2.gif "discrete Y32 spherical harmonic target with varying amount of deformation")

Aiming for target spectrum

![ani#2](https://raw.githubusercontent.com/levincoolxyz/invspec/master/i3_300_t2_abs(Y32(v))_e0.1p0.5-2.gif "discrete Y32 spherical harmonic target with varying amount of deformation")

Testing with bunny

*coming soon*
