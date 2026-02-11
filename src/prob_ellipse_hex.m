function P = prob_ellipse_hex(p, a, b)
%PROB_ELLIPSE_HEX Single-attempt passage probability for an ellipse through a hexagonal aperture.
%
%   P = PROB_ELLIPSE_HEX(p, a, b) computes the single-attempt passage
%   probability of an in-plane ellipse with semi-axes (a,b) through a
%   regular hexagonal opening (centered at the origin) with apothem p.
%
%   The probability is evaluated as the expected admissible centroid area
%   (normalized by the hexagon area) under a uniform orientation model:
%
%       P = (1 / A_H) * \int_0^\pi A(theta) g(theta) dtheta,
%
%   where A(theta) is the admissible area returned by AREA_DTHETA_HEX and
%   g(theta) = 1/pi is the uniform density on [0, pi).
%
%   Inputs:
%     p - hexagon apothem
%     a - ellipse semi-major axis
%     b - ellipse semi-minor axis
%
%   Output:
%     P - single-attempt passage probability
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

A_H = 2*sqrt(3)*p^2;   % hexagon area (apothem p)
g = @(theta) 1/pi;     % uniform orientation density on [0, pi)
integrand = @(theta) g(theta) .* arrayfun(@(t) area_Dtheta_hex(p,a,b,t), theta);
P = (1/A_H) * integral(integrand, 0, pi, 'AbsTol',1e-8,'RelTol',1e-6);

end