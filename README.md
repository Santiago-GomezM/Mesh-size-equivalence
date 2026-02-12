# Probability–based size equivalence for screening apertures

This repository contains the MATLAB code and supplementary material associated
with the paper:

**“A geometric–probabilistic framework for size equivalence in mineral screening”**

This work has been developed within the framework of the AVANTIS project (Grant Agreement No. 101137552). AVANTIS addresses the challenge that Europe hosts numerous unexploited, low-grade vanadium-bearing titanomagnetite deposits located in Finland, Sweden, Greenland, Norway, Poland, and Ukraine. These deposits are characterized by complex, “spiderweb-like” mineral assemblages which make conventional extraction routes economically unfeasible. Without selective blasting, selective fragmentation, and efficient pre-concentration technologies capable of separating Ti-rich ilmenite from V-bearing magnetite, their exploitation remains limited.

The geometric–probabilistic framework and MATLAB implementations provided in this repository contribute to the objectives of AVANTIS by enabling rigorous fragment size characterization and improved interpretation of screening-based measurements. Such quantitative tools are essential for assessing fragmentation quality, mineral liberation, and the downstream efficiency of pre-concentration and beneficiation processes within the project.

## Overview

In screening–based size characterization, particle size is inferred from the
probability of passage through screen openings of nominal size. In industrial
practice, however, screening circuits often combine apertures of different
geometries (square, rectangular, hexagonal, circular), which are not
geometrically equivalent when particles are anisotropic.

This work develops a geometric and probabilistic framework to:
- Model particle passage as a function of particle shape, orientation, and
  aperture geometry,
- Compute passage probabilities for ellipsoidal particles,
- Define probability–preserving size equivalence between different opening
  geometries,
- Correct fragment size distributions obtained from heterogeneous screening
  systems.

Particles are modeled as triaxial ellipsoids with uniformly distributed
orientations, and apertures are represented as convex planar domains. Passage
is formulated as a containment problem based on the admissible center domain
derived through Minkowski erosion.

## Repository structure

- The `src/` directory contains all MATLAB implementations used in the article.
- Each function is documented at a high level in `src/README.md`.

## Reproducibility

All numerical results reported in the paper were generated using the MATLAB
routines provided in this repository. The code is released as supplementary
material to ensure reproducibility and to facilitate further application of
the proposed methodology.

## Citation

If you use this code, please cite the associated article. 
[Zenodo DOI](https://doi.org/10.5281/zenodo.18623011)

