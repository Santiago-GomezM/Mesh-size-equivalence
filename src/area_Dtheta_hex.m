function A = area_Dtheta_hex(p, a, b, theta)
%AREA_DTHETA_HEX Effective admissible area for a hexagonal aperture.
%
%   A = AREA_DTHETA_HEX(p, a, b, theta) returns the area of the set
%   D_theta = H ⊖ E_theta, where H is a regular hexagon centered at the
%   origin with apothem p, and E_theta is an ellipse with semi-axes (a,b)
%   rotated by an angle theta (radians).
%
%   The set D_theta represents the admissible region of particle centroids
%   governing passage through the opening under the geometric–probabilistic
%   framework developed in the accompanying manuscript.
%
%   Inputs:
%     p     - hexagon apothem
%     a, b  - ellipse semi-axes (a: semi-major axis, b: semi-minor axis)
%     theta - ellipse orientation (radians)
%
%   Output:
%     A - area of D_theta (zero if the set is empty)
%
%   Implementation note:
%   The admissible polygon is obtained by successive clipping against the
%   six supporting half-spaces of the eroded hexagon. No early rejection
%   based solely on r_k <= 0 is performed, since feasibility depends on the
%   joint intersection of all constraints.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

    %--- 1) Hexagon outward normals (6 sides) ---
    phis    = pi/6 + (0:5)'*pi/3;
    normals = [cos(phis), sin(phis)];   % 6 x 2

    %--- 2) Rotated ellipse matrix ---
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
    M = R * diag([a^2, b^2]) * R.';     % 2x2

    %--- 3) Support function h_{E_theta}(n_k) ---
    h = zeros(6,1);
    for k = 1:6
        u = normals(k,:).';            % 2x1 column
        h(k) = sqrt(u' * M * u);
    end

    %--- 4) Half-spaces n_k · x <= r_k(theta) ---
    r = p - h;

    %--- 5) Large initial polygon containing all possible D_theta ---
    % Circumradius of the hexagon:
    R_hex = 2*p/sqrt(3);
    L = R_hex + max(a,b);              % safety margin
    P = [ -L -L;
           L -L;
           L  L;
          -L  L ];                     % large square

    %--- 6) Successive clipping by each half-space n_k·x <= r_k ---
    for k = 1:6
        nk = normals(k,:);
        rk = r(k);
        P = clip_polygon_halfspace(P, nk, rk);
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
    % P: Nx2 polygon vertices (CW or CCW order)
    % n: 1x2 outward normal vector
    % r: scalar threshold
    %
    % Returns Q: clipped polygon (possibly empty).

    if isempty(P)
        Q = [];
        return;
    end

    Q = [];
    N = size(P,1);

    for i = 1:N
        A = P(i,:);
        B = P(mod(i,N)+1,:);   % next vertex (cyclic)

        fA = dot(n,A) - r;
        fB = dot(n,B) - r;

        insideA = (fA <= 0);
        insideB = (fB <= 0);

        if insideA && insideB
            Q = [Q; B];

        elseif insideA && ~insideB
            t = fA / (fA - fB);
            I = A + t*(B - A);
            Q = [Q; I];

        elseif ~insideA && insideB
            t = fA / (fA - fB);
            I = A + t*(B - A);
            Q = [Q; I; B];

        else
            % both outside: add nothing
        end
    end
end