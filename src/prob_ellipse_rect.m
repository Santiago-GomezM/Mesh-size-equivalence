function P = prob_ellipse_rect(a,b,rect)
%PROB_ELLIPSE_RECT Single-attempt passage probability for an ellipse through a rectangular aperture.
%
%   P = PROB_ELLIPSE_RECT(a, b, rect) computes the single-attempt passage
%   probability of an in-plane ellipse with semi-axes (a,b) through an
%   axis-aligned rectangular opening defined by rect = [xmin xmax ymin ymax].
%
%   The probability is evaluated as the expected admissible centroid area
%   (normalized by the rectangle area) under a uniform orientation model:
%
%       P = (1 / A_R) * \int_0^\pi A(theta) g(theta) dtheta,
%
%   where A(theta) is the admissible area returned by AREA_DTHETA_RECT and
%   g(theta) = 1/pi is the uniform density on [0, pi).
%
%   Inputs:
%     a    - ellipse semi-major axis
%     b    - ellipse semi-minor axis
%     rect - rectangle bounds [xmin xmax ymin ymax]
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

g = @(theta) 1/pi;            % uniform orientation density on [0, pi)
A_rect = (rect(2)-rect(1)) * (rect(4)-rect(3));
integrand = @(th) g(th) .* arrayfun(@(t) area_Dtheta_rect(rect,a,b,t), th);
P = (1/A_rect) * integral(integrand, 0, pi, 'AbsTol',1e-8,'RelTol',1e-6);

end