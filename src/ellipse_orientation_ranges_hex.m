function [out, int_def] = ellipse_orientation_ranges_hex(p, a, b, C)
%ELLIPSE_ORIENTATION_RANGES_HEX Feasible orientation fraction for an ellipse within a hexagonal aperture.
%
%   [out, int_def] = ELLIPSE_ORIENTATION_RANGES_HEX(p, a, b, C) computes,
%   for each candidate center in C, the fraction of ellipse orientations
%   theta ∈ [0, pi) for which an ellipse of semi-axes (a,b) fits inside a
%   regular hexagon (centered at the origin) with apothem p.
%
%   The feasible set of orientations is returned both as:
%     out     - fraction of admissible orientations in [0, pi)
%               (i.e., |Theta| / pi), one value per center.
%     int_def - cell array of admissible orientation intervals in degrees,
%               each entry containing an Mx2 array with rows [theta_min theta_max]
%               in [0, 180]. An empty array indicates no admissible orientations.
%
%   Inputs:
%     p - hexagon apothem
%     a - ellipse semi-major axis
%     b - ellipse semi-minor axis
%     C - Nx2 array of candidate centers [cx_i, cy_i]
%
%   Output:
%     out     - Nx1 vector with the admissible orientation fraction
%     int_def - Nx1 cell array with admissible orientation intervals (degrees)
%
%   Implementation note:
%   The constraints are derived from the support function of the rotated ellipse
%   along the six outward normals of the hexagon. For each face, the feasible
%   orientations are computed analytically and then intersected across all faces.
%   A dedicated branch handles the nearly circular case a ≈ b.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

    N = size(C,1);

    % Preallocations
    out     = zeros(N,1);
    int_def = cell(N,1);

    % Hexagon outward normals (fixed)
    phis = pi/6 + (0:5)*pi/3;   % 1x6
    ux   = cos(phis);           % 1x6
    uy   = sin(phis);           % 1x6

    % Ellipse parameters
    A2 = a*a;
    B2 = b*b;
    denom = A2 - B2;
    tol = 1e-12;

    % Precompute normal arrays (for broadcasting)
    ux = reshape(ux,1,6);       % 1x6
    uy = reshape(uy,1,6);
    phis = reshape(phis,1,6);

    % Projections d = ux*cx + uy*cy  -> Nx6 (fully vectorized)
    d = C(:,1).*ux + C(:,2).*uy;      % Nx6
    rhs = p - d;                      % Nx6

    % Centers for which no orientation can fit
    invalid = any(rhs <= 0, 2);
    out(invalid) = 0;
    int_def(invalid) = {[]};

    % Nearly circular ellipse (a ≈ b): feasibility independent of theta
    if abs(denom) < 1e-15
        r = sqrt(A2);

        ok = ~invalid & all(rhs >= r - tol, 2);
        out(ok) = 1;
        int_def(ok) = { [0 180] };

        return
    end

    % K for all centers and faces -> Nx6
    K = (2*rhs.^2 - (A2 + B2)) ./ denom;

    %-----------------------------------------------------------
    % Loop over centers; geometric quantities already vectorized
    %-----------------------------------------------------------
    for i = 1:N

        if invalid(i)
            continue
        end

        Ki   = K(i,:);       % 1x6
        rhsi = rhs(i,:);     % 1x6

        % Current feasible interval set (radians), initialized to [0, pi)
        current = [0, pi];

        for k = 1:6

            Kk  = Ki(k);
            phi = phis(k);

            if Kk < -1
                current = [];
                break

            elseif Kk >= 1
                continue

            else
                alpha0 = acos(Kk);
                alpha1 = 0.5*alpha0;
                alpha2 = pi - alpha1;

                th1 = phi - alpha2;
                th2 = phi - alpha1;

                a1 = mod(th1, pi);
                b1 = mod(th2, pi);

                if a1 <= b1
                    allowed = [a1, b1];
                else
                    allowed = [0, b1;
                               a1, pi];
                end

                current = intersect_interval_sets(current, allowed);
                if isempty(current)
                    break
                end
            end
        end

        if isempty(current)
            out(i) = 0;
            int_def{i} = [];
        else
            current = merge_intervals(current);
            L = current(:,2) - current(:,1);
            out(i) = sum(L)/pi;
            int_def{i} = current * 180/pi;
        end
    end
end

%---------------------- AUXILIARY FUNCTIONS ---------------------------%

function S = intersect_interval_sets(A, B)
%INTERSECT_INTERVAL_SETS Intersection of two interval sets on [0, pi).
%
%   S = INTERSECT_INTERVAL_SETS(A, B) returns the intersection of interval
%   sets A and B, where each is an Mx2 (or Nx2) array of intervals [L, R]
%   in radians within [0, pi). The output S is an array of intervals
%   (possibly empty). Overlapping/adjacent intervals are merged.

    if isempty(A) || isempty(B)
        S = [];
        return;
    end

    S = [];
    for i = 1:size(A,1)
        a1 = A(i,1); b1 = A(i,2);
        for j = 1:size(B,1)
            a2 = B(j,1); b2 = B(j,2);
            L = max(a1,a2);
            R = min(b1,b2);
            if R > L
                S = [S; L, R];
            end
        end
    end

    if isempty(S)
        return;
    end

    % Merge potential overlaps
    S = merge_intervals(S);
end

function I = merge_intervals(I)
%MERGE_INTERVALS Merge overlapping or adjacent intervals.
%
%   I = MERGE_INTERVALS(I) merges overlapping or numerically adjacent
%   intervals (within a small tolerance) in an Mx2 array of intervals [L, R].

    if isempty(I)
        return;
    end
    I = sortrows(I,1);
    res = I(1,:);
    for k = 2:size(I,1)
        a = I(k,1); b = I(k,2);
        if a <= res(end,2) + 1e-12  % overlap or touch (with tolerance)
            res(end,2) = max(res(end,2), b);
        else
            res = [res; a, b];
        end
    end
    I = res;
end