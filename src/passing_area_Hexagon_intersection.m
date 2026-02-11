function [X, Y, percepass] = passing_area_Hexagon_intersection(p, a, b, nTheta, nPhi)
%PASSING_AREA_HEXAGON_INTERSECTION Centroid region valid for all ellipse orientations in a hexagonal aperture.
%
%   [X, Y, percepass] = PASSING_AREA_HEXAGON_INTERSECTION(p, a, b, nTheta, nPhi)
%   computes the planar region of ellipse centers for which an ellipse with
%   semi-axes (a,b) fits inside a regular hexagon (centered at the origin,
%   apothem p) for all in-plane orientations theta.
%
%   The region is constructed in polar form by evaluating, for each ray
%   direction phi, the maximum admissible radial distance r(phi) obtained
%   from the most restrictive combination of:
%     (i) hexagon side constraints (six outward normals), and
%     (ii) the worst-case ellipse orientation theta.
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

% Calcula la región donde la elipse E(a,b) cabe en el hexágono
% para todas las orientaciones theta ∈ [0, 2π].

if nargin < 4, nTheta = 500; end
if nargin < 5, nPhi = 500; end

alpha = pi/6 + (0:5)*pi/3;      % Normales del hexágono (6)
theta = linspace(0, 2*pi, nTheta);

% Precompute support function h_k(theta)
ang = alpha' - theta;           % 6 x nTheta
Delta = sqrt((a*cos(ang)).^2 + (b*sin(ang)).^2);   % h_k(theta)

phi = linspace(0, 2*pi, nPhi);
r = zeros(1, nPhi);

for j = 1:nPhi
    denom = abs(cos(phi(j) - alpha));   % 1 x 6
    denom(denom < 1e-12) = NaN;         % evitar 0
    denomMat = repmat(denom', 1, nTheta);

    aux = p - Delta;       % p - h_k(theta)
    aux(aux < 0) = 0;      % si ya no cabe, radio = 0

    limits = aux ./ denomMat;   % 6 x nTheta

    % Intersection (worst-case over orientations)
    % First: minimum over constraints for each theta (most restrictive side)
    minLimits = min(limits, [], 1);     % 1 x nTheta

    % Then: minimum over theta (worst orientation)
    r(j) = min(minLimits);
end

% Convert to cartesian coordinates
X = r .* cos(phi);
Y = r .* sin(phi);

areaPermissive = polyarea(X, Y);
hexArea = 2*sqrt(3)*p^2;
percepass = areaPermissive / hexArea;

end