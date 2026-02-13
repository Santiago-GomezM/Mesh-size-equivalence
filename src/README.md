# MATLAB source code — Probability of passage and size equivalence

This directory contains the MATLAB source code used in the article  
“A geometric–probabilistic framework for size equivalence in mineral screening”.

The routines implement the geometric and probabilistic framework developed
to compute particle passage probabilities through screen openings of different
geometries and to establish probability–preserving size equivalence between
apertures.

Particles are modeled as ellipsoids with random orientations, and openings are
represented as convex planar domains (hexagonal, rectangular, square, and
circular).

---

## 1. Admissible center domains for projected ellipses

These functions compute the area of the admissible center domain \(D_\theta\)
associated with a projected ellipse of fixed in–plane orientation.

### Projection from a 3D ellipsoid
- **`area_Dabg_hex.m`**  
  Computes the area of D<sub>αβγ</sub> for an ellipsoid passing through a
  regular hexagonal opening for a given set of Euler angles.

- **`area_Dabg_rect.m`**  
  Computes the area of D<sub>αβγ</sub> for an ellipsoid passing through a
  hexagonal hexagonal opening for a given set of Euler angles.

### Fixed in–plane orientation
- **`area_Dtheta_hex.m`**  
  Area of D<sub>θ</sub> for a rotated ellipse inside a regular hexagon.

- **`area_Dtheta_rect.m`**  
  Area of D<sub>θ</sub> for a rotated ellipse inside a rectangle.

- **`area_Dtheta_circ.m`**  
  Area of D<sub>θ</sub> for a circular opening (orientation–independent).

---

## 2. Restrictive and permissive regions

These routines compute the boundaries of the intersection and union of
orientation–dependent admissible regions, corresponding respectively to passage
for all orientations and for at least one orientation.

### Hexagonal openings
- **`passing_area_Hexagon_intersection.m`**  
  Restrictive region |D<sup>∩</sup>| for hexagonal openings.

- **`passing_area_Hexagon_union.m`**  
  Permissive region |D/∪</sup>| for hexagonal openings.

### Rectangular openings
- **`passing_area_Rectangle_intersection.m`**  
  Restrictive region |D<sup>∩</sup>| for rectangular openings.

- **`passing_area_Rectangle_union.m`**  
  Permissive region |D<sup>∪</sup>| for rectangular openings.

---

## 3. Feasible orientation ranges

These functions compute the set of in–plane orientations for which an ellipse
centered at a given position can pass through an opening.

- **`ellipse_orientation_ranges_hex.m`**  
  Feasible orientation intervals for hexagonal openings.

- **`ellipse_orientation_ranges_rect.m`**  
  Feasible orientation intervals for rectangular openings.

---

## 4. Probability of passage for projected ellipses (2D)

- **`prob_ellipse_hex.m`**  
  Single–attempt probability of passage through a hexagonal opening.

- **`prob_ellipse_rect.m`**  
  Single–attempt probability of passage through a rectangular opening.

---

## 5. Probability of passage for ellipsoids (3D orientations)

These routines average over uniformly distributed orientations in \(SO(3)\).

- **`prob_ellipsoid_hex.m`**  
  Probability of passage of an ellipsoid through a hexagonal opening.

- **`prob_ellipsoid_rect.m`**  
  Probability of passage of an ellipsoid through a rectangular opening.

---

## 6. Size equivalence between openings

- **`equivalent_opening_hex.m`**  
  Computes probability–preserving size ratios between hexagonal and square
  openings.

- **`equivalent_opening_rect.m`**  
  Computes probability–preserving size ratios between rectangular and square
  openings.

---

## 7. Maximum fragment size based on equivalent square openings

- **`equivalent_square.m`**  
  Determines the maximum fragment size based on probabilistic passage criteria
  through an equivalent square opening.

---

## Notes
- All routines are written in MATLAB and are self–contained.
- No external toolboxes are required.
- Angles are expressed in radians unless otherwise stated.
