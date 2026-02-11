function [X, Y, rI] = passing_area_Rectangle_intersection(a,b,rect,nTheta,nPsi)
%PASSING_AREA_RECTANGLE_INTERSECTION Centroid region valid for all ellipse orientations in a rectangular aperture.
%
%   [X, Y, rI] = PASSING_AREA_RECTANGLE_INTERSECTION(a, b, rect, nTheta, nPsi)
%   computes the planar region of ellipse centers for which an ellipse with
%   semi-axes (a,b) fits inside an axis-aligned rectangle rect = [xmin xmax ymin ymax]
%   for all in-plane orientations theta (intersection over theta).
%
%   The region is constructed in polar form by evaluating, for each ray
%   direction psi, a radial limit rI(psi) corresponding to the worst-case
%   ellipse orientation. The result is then clipped by the rectangle boundary.
%
%   Outputs:
%     X, Y - boundary coordinates of the admissible region (cartesian)
%     rI   - admissible radial distances rI(psi) for each polar direction psi
%
%   Inputs:
%     a, b   - ellipse semi-axes
%     rect   - rectangle bounds [xmin xmax ymin ymax]
%     nTheta - number of sampled orientations theta in [0, pi] (default: 400)
%     nPsi   - number of sampled polar directions psi in [-pi, pi] (default: 720)
%
%   Notes:
%   - The returned boundary is an approximation whose accuracy improves with
%     nTheta and nPsi.
%   - Orientation is sampled in [0, pi] due to pi-periodicity of the ellipse.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

    if nargin < 4, nTheta = 400; end
    if nargin < 5, nPsi   = 720; end

    % Polar directions and ellipse orientations
    psis   = linspace(-pi, pi, nPsi);
    thetas = linspace(0,  pi, nTheta);

    % Rectangle geometry (outward normals and offsets)
    xmin=rect(1); xmax=rect(2); ymin=rect(3); ymax=rect(4);
    normals = [1 0; -1 0; 0 1; 0 -1];
    d       = [xmax; -xmin; ymax; -ymin];

    A2 = a*a; B2 = b*b;

    rI = zeros(size(psis));

    for j = 1:nPsi
        psi = psis(j);

        % Maximum radius inside the rectangle along direction psi
        rmax = boundary_rect_polar(rect, psi);

        % If the ray does not intersect the rectangle, radius is zero
        if rmax <= 0
            rI(j) = 0;
            continue;
        end

        % rI(psi) = min_theta r_theta(psi) (worst-case over orientations)
        rmin = inf;

        for t = 1:nTheta
            th = thetas(t);

            % Ellipse support along the 4 rectangle normals
            h = zeros(4,1);
            for k = 1:4
                phi_k = atan2(normals(k,2), normals(k,1));
                h(k) = sqrt(A2*cos(phi_k - th)^2 + B2*sin(phi_k - th)^2);
            end

            % Radial limits imposed by each side (when the side constrains direction psi)
            r_k = inf(4,1);
            for k = 1:4
                phi_k = atan2(normals(k,2), normals(k,1));
                cp = cos(psi - phi_k);
                if cp > 1e-12
                    r_k(k) = (d(k) - h(k)) / cp;
                end
            end

            % If any side yields a non-positive bound, no admissible radius for this theta
            if any(r_k <= 0)
                r_theta = 0;
            else
                r_theta = min(r_k);
            end

            rmin = min(rmin, r_theta);
        end

        % Clip by the rectangle boundary
        rI(j) = max(0, min(rmin, rmax));
    end

    X = rI .* cos(psis);
    Y = rI .* sin(psis);
end

%====================================================================
% Maximum radius up to the rectangle boundary along direction psi
%====================================================================
function rmax = boundary_rect_polar(rect, psi)
    xmin=rect(1); xmax=rect(2); ymin=rect(3); ymax=rect(4);
    cpsi = cos(psi); spsi = sin(psi);

    candidates = [];

    % Intersections with x = xmin and x = xmax
    if abs(cpsi) > 1e-12
        r1 = xmin / cpsi;
        r2 = xmax / cpsi;
        candidates = [candidates, r1, r2];
    end

    % Intersections with y = ymin and y = ymax
    if abs(spsi) > 1e-12
        r3 = ymin / spsi;
        r4 = ymax / spsi;
        candidates = [candidates, r3, r4];
    end

    % Keep positive radii whose intersection point lies within the rectangle
    r_valid = [];
    for r = candidates
        if r > 0
            x = r*cpsi; y = r*spsi;
            if x >= xmin-1e-9 && x <= xmax+1e-9 && ...
               y >= ymin-1e-9 && y <= ymax+1e-9
                r_valid(end+1) = r; 
            end
        end
    end

    if isempty(r_valid)
        rmax = 0;
    else
        rmax = min(r_valid);
    end
end