function A = area_Dtheta_rect(rect, a, b, theta)
%AREA_DTHETA_RECT Effective admissible area for a rectangular aperture.
%
%   A = AREA_DTHETA_RECT(rect, a, b, theta) returns the area of the set
%   D_theta = R ⊖ E_theta, where R is an axis-aligned rectangle and E_theta
%   is an ellipse with semi-axes (a,b) rotated by an angle theta (radians).
%
%   The set D_theta represents the admissible region of particle centroids
%   governing passage through the opening under the geometric–probabilistic
%   framework developed in the accompanying manuscript.
%
%   Inputs:
%     rect  - rectangle bounds [xmin xmax ymin ymax]
%     a, b  - ellipse semi-axes (typically a >= b)
%     theta - ellipse orientation (radians)
%
%   Output:
%     A - area of D_theta (zero if the set is empty)
%
%   Implementation note:
%   The admissible polygon is obtained by successive clipping against the
%   four supporting half-spaces of the eroded rectangle.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

    xmin = rect(1); xmax = rect(2);
    ymin = rect(3); ymax = rect(4);

    %--- 1) Rectangle outward normals (sides) ---
    normals = [ 1 0;   % x <= xmax
               -1 0;   % x >= xmin
                0 1;   % y <= ymax
                0 -1]; % y >= ymin

    d = [xmax; -xmin; ymax; -ymin];

    %--- 2) Rotated ellipse matrix ---
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    M = R * diag([a*a, b*b]) * R.';

    %--- 3) Support function h_{E_theta}(n_k) ---
    h = zeros(4,1);
    for k = 1:4
        u = normals(k,:)';
        h(k) = sqrt(u' * M * u);
    end

    %--- 4) Half-spaces n_k · x <= r_k(theta) ---
    r = d - h;

    %--- 5) Large initial polygon containing all possible D_theta ---
    %     Use a sufficiently large square that contains D_theta if non-empty.
    L = max([abs(xmin), abs(xmax), abs(ymin), abs(ymax)]) + max(a,b);
    P = [ -L -L;
           L -L;
           L  L;
          -L  L ];

    %--- 6) Successive clipping by the 4 half-spaces ---
    for k = 1:4
        nk = normals(k,:);
        rk = r(k);
        P  = clip_polygon_halfspace(P, nk, rk);
        if isempty(P)
            A = 0;
            return;
        end
    end

    %--- 7) Area of the resulting polygon ---
    A = polyarea(P(:,1), P(:,2));
end


%==================================================================
% Clip polygon P by the half-space { x : n·x <= r }
%==================================================================
function Q = clip_polygon_halfspace(P, n, r)
    if isempty(P)
        Q = [];
        return;
    end

    Q = [];
    N = size(P,1);

    for i = 1:N
        A = P(i,:);
        B = P(mod(i,N)+1,:);

        fA = dot(n,A) - r;
        fB = dot(n,B) - r;

        insideA = (fA <= 0);
        insideB = (fB <= 0);

        if insideA && insideB
            Q = [Q; B];
        elseif insideA && ~insideB
            t = fA / (fA - fB);
            I = A + t*(B-A);
            Q = [Q; I];
        elseif ~insideA && insideB
            t = fA / (fA - fB);
            I = A + t*(B-A);
            Q = [Q; I; B];
        end
    end
end