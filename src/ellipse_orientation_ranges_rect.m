function [out, int_def] = ellipse_orientation_ranges_rect(rect,a,b,C)
%ELLIPSE_ORIENTATION_RANGES_RECT Feasible orientation fraction for an ellipse within a rectangular aperture.
%
%   [out, int_def] = ELLIPSE_ORIENTATION_RANGES_RECT(rect, a, b, C) computes,
%   for each candidate center in C, the fraction of ellipse orientations
%   theta ∈ [0, pi) for which an ellipse of semi-axes (a,b) fits inside an
%   axis-aligned rectangle defined by rect = [xmin xmax ymin ymax].
%
%   The feasible set of orientations is returned both as:
%     out     - fraction of admissible orientations in [0, pi)
%               (i.e., |Theta| / pi), one value per center.
%     int_def - cell array of admissible orientation intervals in degrees,
%               each entry containing an Mx2 array with rows [theta_min theta_max]
%               in [0, 180]. An empty array indicates no admissible orientations.
%
%   Inputs:
%     rect - rectangle bounds [xmin xmax ymin ymax]
%     a    - ellipse semi-major axis
%     b    - ellipse semi-minor axis
%     C    - Nx2 array of candidate centers [cx_i, cy_i]
%
%   Output:
%     out     - Nx1 vector with the admissible orientation fraction
%     int_def - Nx1 cell array with admissible orientation intervals (degrees)
%
%   Implementation note:
%   The constraints are derived from the support function of the rotated ellipse
%   along the four outward normals of the rectangle. For each side, the feasible
%   orientations are computed analytically and then intersected across all sides.
%   A dedicated branch handles the nearly circular case a ≈ b.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

    N = size(C,1);

    %--------------------------------------------------
    % Preallocation
    %--------------------------------------------------
    out     = zeros(N,1);
    int_def = cell(N,1);

    %--------------------------------------------------
    % Fixed rectangle geometry
    %--------------------------------------------------
    normals = [ 1  0;    % x <= xmax
               -1  0;    % x >= xmin
                0  1;    % y <= ymax
                0 -1 ];  % y >= ymin

    d = [rect(2); -rect(1); rect(4); -rect(3)];   % offsets for each side

    phis = [0; pi; pi/2; -pi/2];                  % normal directions (radians)

    %--------------------------------------------------
    % Ellipse parameters
    %--------------------------------------------------
    A2 = a*a;
    B2 = b*b;
    denom = A2 - B2;
    tol = 1e-12;

    %--------------------------------------------------
    % Vectorized geometric part
    %--------------------------------------------------
    % r(i,k) = available margin for center i with respect to side k
    r = d.' - C*normals.';        % Nx4

    % Centers for which no orientation can fit
    invalid = any(r <= 0, 2);
    out(invalid) = 0;
    int_def(invalid) = {[]};

    %--------------------------------------------------
    % Nearly circular case (a ≈ b): feasibility independent of theta
    %--------------------------------------------------
    if abs(denom) < 1e-15
        h = sqrt(A2);
        ok = ~invalid & all(r >= h - tol, 2);
        out(ok) = 1;
        int_def(ok) = { [0 180] };
        return
    end

    %--------------------------------------------------
    % K coefficients for all centers and sides
    %--------------------------------------------------
    K = (2*r.^2 - (A2 + B2)) ./ denom;   % Nx4

    %--------------------------------------------------
    % Interval logic loop over centers
    %--------------------------------------------------
    for i = 1:N

        if invalid(i)
            continue
        end

        Ki = K(i,:);    % 1x4

        % Initial interval set: full [0, pi)
        current = [0 pi];

        for k = 1:4

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

                % theta ∈ [phi - alpha2, phi - alpha1]
                th1 = phi - alpha2;
                th2 = phi - alpha1;

                a1 = mod(th1, pi);
                b1 = mod(th2, pi);

                if a1 <= b1
                    allowed = [a1 b1];
                else
                    allowed = [0 b1;
                               a1 pi];
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
            out(i) = sum(L) / pi;
            int_def{i} = current * 180/pi;
        end
    end
end

%====================================================================
% Interval set utilities on [0, pi)
%====================================================================
function S = intersect_interval_sets(A,B)
%INTERSECT_INTERVAL_SETS Intersection of two interval sets on [0, pi).
%
%   S = INTERSECT_INTERVAL_SETS(A, B) returns the intersection of interval
%   sets A and B, where each is an Mx2 (or Nx2) array of intervals [L, R]
%   in radians within [0, pi). The output S is an array of intervals
%   (possibly empty). Overlapping/adjacent intervals are merged.

    if isempty(A) || isempty(B)
        S = [];
        return
    end
    S = [];
    for i = 1:size(A,1)
        for j = 1:size(B,1)
            L = max(A(i,1), B(j,1));
            R = min(A(i,2), B(j,2));
            if R > L
                S = [S; L R]; %#ok<AGROW>
            end
        end
    end
    S = merge_intervals(S);
end

function I = merge_intervals(I)
%MERGE_INTERVALS Merge overlapping or adjacent intervals.
%
%   I = MERGE_INTERVALS(I) merges overlapping or numerically adjacent
%   intervals (within a small tolerance) in an Mx2 array of intervals [L, R].

    if isempty(I), return; end
    I = sortrows(I,1);
    res = I(1,:);
    for k = 2:size(I,1)
        if I(k,1) <= res(end,2) + 1e-12
            res(end,2) = max(res(end,2), I(k,2));
        else
            res = [res; I(k,:)]; %#ok<AGROW>
        end
    end
    I = res;
end