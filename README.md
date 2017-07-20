# 2D Discrete Inverse Spectral Problem

### Goal 
To calculate an approximate solution to the classic (Laplace-Beltrami) inverse spectral problem for discrete (genus 0) surfaces.

### Content
Included is a suite of MATLAB codes implementing the naive direct gradient descent approach.

test_script.m is the top-level script that generates results.

### People
Project envisioned, advised, and supervised by **Prof. Etienne Vouga** and **Prof. Keenan Crane**

Some codes here (on mesh optimization and a demo of spherical harmonics) are not mine.

### Testing Procedures
1. use conformalized mean curvature flow (cMCF) to get a spherical mesh with a target set of conformal factors from a target mesh

2. use BFGS descent search for some conformal factors that achieve a spectrum similar to the one desired

3. embed the metric to a resulting mesh from the sphere by optimizing edge lengths obtained from the factors

4. compare with the target mesh, target spectrum, and cMCF conformal factors

### Results

![ani#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i3_300_t2_abs(Y33(v))_e0.1-1p0.5.gif "discrete Y33 spherical harmonic target with varying percent of eigenvalues used")
![ani#2](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i3_300_t2_abs(Y32(v))_e0.1p0.5-2.gif "discrete Y32 spherical harmonic target with varying amount of deformation")

(cheating) optimize for cMCF spectrum instead
![ani#3](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i2_300_t2_abs(Y32(v))_e0.1p0.5-2.gif "discrete Y32 spherical harmonic target with varying amount of deformation")

**Tests with Spot the cow**

smoothed cow without bi-laplacian regularization
![spot#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i4_mcf_t4_cow03_e1p0r0.png "smoothed Spot as target without regularization")

original spot with regularization

![spot#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i4_mcf_t3_spot1k_e1p0r0.1.png "Spot as target with regularization")

**Tests with bunny**

![bun#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i4_mcf_t3_bunny327_e0.95p0r0.05.png "classic bunny as target with regularization")

*recursive run*
![bun#2](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i3_bun3_t3_bunny327_e0.5p0r0.01.png "classic bunny as target with regularization fitted recursively")

*finer mesh recursive run*
![bun#3](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i3_bun2_t3_bunny602_e0.95p0r0.1.png "classic bunny as target with regularization fitted recursively")

*Before and After* (with minor smoothing)

![bbfe](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/before.png)![baft](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/after.png)

*even finer mesh recursive run*
![bun#3](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i3_bun2_t3_bunny1043_e0.5p0r0.05.png "classic bunny as target with regularization fitted recursively")

**"Mesh-free" Spherical Harmonic Basis Solution**

From now on we have *number of eigenvalues used = number of free SH basis function coefficient = n, LB operator expanded in 961 SH basis functions*

PL spectrum as target: n = 36
![blob#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i2_540_t3_blob18k_a36e36L30.png "blob mesh PL spectrum as target")
n = 49
![blob#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i2_540_t3_blob18k_a49e49L30.png "blob mesh PL spectrum as target")
n = 64
![blob#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i2_540_t3_blob18k_a64e64L30.png "blob mesh PL spectrum as target")

SH spectrum as target (cheating): n = 49
![blob#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/doc/page/i2_540_t3_blob18k_a49e49L30s.png "blob mesh PL spectrum as target")

Results were adjusted up to SO(3) to mod out the rigid rotation ambiguity
### Ways to go

1. (ongoing) without prior knowledge of the target mesh, we will have to start from a uniform (coarse) spherical mesh and develop a suitable adaptive refinement scheme

2. (banging my head) why does high frequency data matter in the FEM/hat function basis? the current way involves optimizing for them and then penalize for its noisyness via regularization, which seems very silly...

2. in practice the inverse problem would not be about the Laplace-Beltrami operator (need to consider bending energy of thin shell etc.)

2. in practice higher frequencies will most definitly be prohibitively noisy

3. can we guess the topology beforehand? would higher genus surfaces work in similar fashion despite planarity/hyperbolicity? (e.g. there are known non-trivial isospectral hyperbolic (g>5) surfaces...)
(Yes. Reuter, Wolter, Peinecke 2006 ~ first 500 eigenvalues)