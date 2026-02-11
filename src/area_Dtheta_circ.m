function A = area_Dtheta_circ(R, a)

%   AREA_DTHETA_CIRC Effective admissible area for a circular aperture.
%
%   A = AREA_DTHETA_CIRC(R, a) computes the admissible area D_theta
%   corresponding to the Minkowski difference B_R ⊖ E_theta, where B_R
%   is a circular aperture of radius R and E_theta is the projected
%   ellipse of semi-major axis a. The result is independent of the
%   orientation angle theta and of the minor semi-axis.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

    if R <= a
        A = 0;
    else
        A = pi * (R - a)^2;
    end
end