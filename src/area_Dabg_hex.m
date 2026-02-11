function A = area_Dabg_hex(p, a3, b3, c3, alpha, beta, gamma)
%AREA_DABG_HEX Effective admissible area for a hexagonal aperture (3D orientation).
%
%   A = AREA_DABG_HEX(p, a3, b3, c3, alpha, beta, gamma) computes the
%   admissible area D_theta for the passage of an ellipsoidal particle
%   through a regular hexagonal opening of apothem p.
%
%   The particle is modeled as a triaxial ellipsoid with semi-axes
%   (a3, b3, c3) and orientation defined by the Euler angles
%   (alpha, beta, gamma). The function:
%       1) constructs the quadratic form of the rotated ellipsoid,
%       2) derives the projected ellipse in the screen plane,
%       3) evaluates the admissible area via AREA_DTHETA_HEX.
%
%   The resulting area represents the centroid domain associated with
%   particle passage under the geometric–probabilistic framework
%   developed in the accompanying manuscript.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

Dinv = diag([1/a3^2, 1/b3^2, 1/c3^2]);

Rz1 = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
Ry  = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rz2 = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];

R = Rz1 * Ry * Rz2;

Q = R * Dinv * R.';

Q2  = Q(1:2,1:2);
q   = Q(1:2,3);
Q33 = Q(3,3);

M = Q2 - (q*q.')/Q33;

% Eigenvalues of symmetric 2x2 matrix M (closed-form solution)
m11 = M(1,1); m12 = M(1,2); m22 = M(2,2);
tr  = m11 + m22;
disc = sqrt((m11 - m22)^2 + 4*m12^2);
lam1 = 0.5*(tr - disc);
lam2 = 0.5*(tr + disc);

a_proj = 1/sqrt(lam1);
b_proj = 1/sqrt(lam2);

theta = 0.5*atan2(2*m12, (m11 - m22));   % in (-pi/2, pi/2]

A = area_Dtheta_hex(p, a_proj, b_proj, theta);

end