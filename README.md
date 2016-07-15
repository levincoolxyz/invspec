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

![ani#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/i3_300_t2_abs(Y33(v))_e0.1-1p0.5.gif "discrete Y33 spherical harmonic target with varying percent of eigenvalues used")
![ani#2](https://raw.githubusercontent.com/levincoolxyz/invspec/master/i3_300_t2_abs(Y32(v))_e0.1p0.5-2.gif "discrete Y32 spherical harmonic target with varying amount of deformation")

(cheating) optimize for cMCF spectrum instead
![ani#3](https://raw.githubusercontent.com/levincoolxyz/invspec/master/i2_300_t2_abs(Y32(v))_e0.1p0.5-2.gif "discrete Y32 spherical harmonic target with varying amount of deformation")

**Tests with Spot the cow**

smoothed cow without bi-laplacian regularization
![spot#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/spot/cow/i4_mcf_t4_cow03_e1p0r0.png "smoothed Spot as target without regularization")

original spot with regularization

![spot#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/spot/i4_mcf_t3_spot1k_e1p0r0.1.png "Spot as target with regularization")

**Tests with bunny**

![bun#1](https://raw.githubusercontent.com/levincoolxyz/invspec/master/bunny/i4_mcf_t3_bunny327_e0.95p0r0.05.png "classic bunny as target with regularization")

*recursive run*
![bun#2](https://raw.githubusercontent.com/levincoolxyz/invspec/master/bunny/recursive327/i3_bun3_t3_bunny327_e0.5p0r0.01.png "classic bunny as target with regularization fitted recursively")

*finer mesh recursive run*
![bun#3](https://raw.githubusercontent.com/levincoolxyz/invspec/master/bunny/recursive602/i3_bun2_t3_bunny602_e0.95p0r0.1.png "classic bunny as target with regularization fitted recursively")

*Before and After* (with minor smoothing)

![bbfe](https://raw.githubusercontent.com/levincoolxyz/invspec/master/bunny/recursive602/before.png)![baft](https://raw.githubusercontent.com/levincoolxyz/invspec/master/bunny/recursive602/after.png)

*even finer mesh recursive run*
![bun#3](https://raw.githubusercontent.com/levincoolxyz/invspec/master/bunny/recursive1k/i3_bun2_t3_bunny1k_e0.5p0r0.05.png "classic bunny as target with regularization fitted recursively")

*more coming soon*

### Ways to go

1. (ongoing) without prior knowledge of the target mesh, we will have to start from a uniform (coarse) spherical mesh and develop a suitable adaptive refinement scheme

2. (banging my head) why does high frequency data matter in the FEM/hat function basis? the current way involves optimizing for them and then penalize for its noisyness via regularization, which seems very silly...

2. (ongoing) one way to [maybe] resolve the above problem is to use the mesh independent "Fourier" spherical harmonic basis, but will it blend?

2. in practice the inverse problem would not be about the Laplace-Beltrami operator (need to consider bending energy of thin shell etc.)

2. in practice higher frequencies will most definitly be prohibitively noisy

3. can we guess the topology beforehand? would higher genus surfaces work in similar fashion despite planarity/hyperbolicity? (e.g. there are known non-trivial isospectral hyperbolic (g>5) surfaces...)