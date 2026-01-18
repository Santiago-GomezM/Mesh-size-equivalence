# Probability–based size equivalence for screening apertures

This repository contains the MATLAB code and supplementary material associated
with the article:

**“[Full title of the paper]”**

## Overview

In screening–based size characterization, particle size is inferred from the
probability of passage through screen openings of nominal size. In industrial
practice, however, screening circuits often combine apertures of different
geometries (square, rectangular, hexagonal, circular), which are not
geometrically equivalent when particles are anisotropic.

This work develops a geometric and probabilistic framework to:
- model particle passage as a function of particle shape, orientation, and
  aperture geometry,
- compute passage probabilities for ellipsoidal particles,
- define probability–preserving size equivalence between different opening
  geometries,
- correct fragment size distributions obtained from heterogeneous screening
  systems.

Particles are modeled as triaxial ellipsoids with uniformly distributed
orientations, and apertures are represented as convex planar domains. Passage
is formulated as a containment problem based on the admissible center domain
derived through Minkowski erosion.

## Repository structure

.
├── src/ # Core MATLAB routines
├── README.md # This file
└── LICENSE

css
Copiar código

- The `src/` directory contains all MATLAB implementations used in the article.
- Each function is documented at a high level in `src/README.md`.

## Reproducibility

All numerical results reported in the paper were generated using the MATLAB
routines provided in this repository. The code is released as supplementary
material to ensure reproducibility and to facilitate further application of
the proposed methodology.

## Citation

If you use this code, please cite the associated article. A DOI-linked archive
(e.g. Zenodo) may be provided for long-term preservation.
