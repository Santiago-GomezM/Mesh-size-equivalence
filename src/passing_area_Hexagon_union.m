function [X, Y, percepass] = passing_area_Hexagon_union(p, a, b, nTheta, nPhi)
%PASSING_AREA_HEXAGON_UNION Centroid region valid for at least one ellipse orientation in a hexagonal aperture.
%
%   [X, Y, percepass] = PASSING_AREA_HEXAGON_UNION(p, a, b, nTheta, nPhi)
%   computes the planar region of ellipse centers for which an ellipse with
%   semi-axes (a,b) fits inside a regular hexagon (centered at the origin,
%   apothem p) for at least one in-plane orientation theta.
%
%   In contrast with PASSING_AREA_HEXAGON_INTERSECTION (worst-case over
%   orientations), this function performs a union over orientations:
%     (i) for each orientation theta, the admissible radial limit is
%         determined by the most restrictive hexagon side,
%     (ii) for each ray direction phi, the maximum radius over theta is taken.
%
%   Outputs:
%     X, Y      - boundary coordinates of the admissible region (cartesian)
%     percepass - admissible area normalized by the hexagon area
%
%   Inputs:
%     p      - hexagon apothem
%     a, b   - ellipse semi-axes
%     nTheta - number of sampled orientations theta in [0, 2*pi] (default: 500)
%     nPhi   - number of sampled polar directions phi in [0, 2*pi] (default: 500)
%
%   Notes:
%   - The support function h_k(theta) is precomputed for the six hexagon
%     outward normals and all sampled orientations.
%   - The returned boundary is an approximation whose accuracy improves with
%     nTheta and nPhi.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

if nargin < 4, nTheta = 500; end
if nargin < 5, nPhi   = 500; end

alpha = pi/6 + (0:5)*pi/3;          % 1×6 hexagon outward normals
theta = linspace(0, 2*pi, nTheta);  % 1×nTheta ellipse orientations

% ---- Precompute support function Delta (6 × nTheta) ----
ang   = alpha.' - theta;            % 6×nTheta
Delta = sqrt((a*cos(ang)).^2 + (b*sin(ang)).^2);

% ---- Vectorized computation for all ray directions phi ----
phi = linspace(0, 2*pi, nPhi);      % 1×nPhi

% Absolute cosine term as in original formulation
denom = abs(cos(phi - alpha.'));    % 6×nPhi
denom(denom < 1e-9) = NaN;          % avoid near-zero division

% Margin between hexagon boundary and ellipse support
aux = p - Delta;                    % 6×nTheta
aux(aux < 0) = 0;                   % negative margins not allowed

% Reshape for implicit 3D expansion:
%   aux3 : 6×nTheta×1
%   den3 : 6×1×nPhi
aux3 = reshape(aux,   6, nTheta, 1);
den3 = reshape(denom, 6, 1,      nPhi);

% limits: 6 × nTheta × nPhi
limits = aux3 ./ den3;

% Minimum over hexagon sides → nTheta×nPhi
minAcrossSides = squeeze(min(limits, [], 1));  % nTheta×nPhi

% Maximum over ellipse orientations → 1×nPhi (union over theta)
r = max(minAcrossSides, [], 1);                % 1×nPhi

% ---- Convert to Cartesian coordinates ----
X = r .* cos(phi);
Y = r .* sin(phi);

areaPermissive = polyarea(X, Y);
hexArea        = 2*sqrt(3)*p^2;
percepass      = areaPermissive / hexArea;

end
