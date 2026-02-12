function [X,Y,rU] = passing_area_Rectangle_union(a,b,rect,nTheta,nPsi)
%PASSING_AREA_RECTANGLE_UNION Region valid for at least one ellipse orientation in a rectangular aperture.
%
%   [X, Y, rU] = PASSING_AREA_RECTANGLE_UNION(a, b, rect, nTheta, nPsi)
%   computes the planar region of ellipse centers for which an ellipse with
%   semi-axes (a,b) fits inside an axis-aligned rectangle
%   rect = [xmin xmax ymin ymax] for at least one in-plane orientation theta
%   (union over theta).
%
%   The region is constructed in polar form by evaluating, for each ray
%   direction psi, the maximum admissible radial distance rU(psi) over the
%   sampled orientations theta. The result is then clipped by the rectangle
%   boundary along the same direction.
%
%   Outputs:
%     X, Y - boundary coordinates of the admissible region (cartesian)
%     rU   - admissible radial distances rU(psi) for each polar direction psi
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

if nargin<4, nTheta=400; end
if nargin<5, nPsi=720; end

phis = linspace(-pi,pi,nPsi);
thetas = linspace(0,pi,nTheta);

xmin=rect(1); xmax=rect(2); ymin=rect(3); ymax=rect(4);

% Rectangle geometry (outward normals and offsets)
normals = [1 0; -1 0; 0 1; 0 -1];
d = [xmax; -xmin; ymax; -ymin];

A2=a*a; B2=b*b;

rU = zeros(size(phis));

for j=1:nPsi
    psi=phis(j);

    % Maximum radius inside the rectangle along direction psi
    rmax=boundary_rect_polar(rect,psi);

    % Best admissible radius over orientations (union over theta)
    best_r=0;

    for t=1:nTheta
        th=thetas(t);

        % Ellipse support along the 4 rectangle normals
        h = zeros(4,1);
        for k=1:4
            phi_k = atan2(normals(k,2), normals(k,1));
            h(k)=sqrt(A2*cos(phi_k-th)^2 + B2*sin(phi_k-th)^2);
        end

        % Radial limits for direction psi (side constraints)
        r_k = (d - h) ./ max(cos(psi - atan2(normals(:,2), normals(:,1))), 1e-12);

        r_k(r_k<=0) = 0;
        r_theta = min(r_k);

        best_r = max(best_r, r_theta);
    end

    % Clip by the rectangle boundary
    rU(j) = min(best_r, rmax);
end

X = rU .* cos(phis);
Y = rU .* sin(phis);
end

%====================================================================
% Maximum radius up to the rectangle boundary along direction psi
%====================================================================
function rmax = boundary_rect_polar(rect,psi)
xmin=rect(1); xmax=rect(2); ymin=rect(3); ymax=rect(4);
cpsi=cos(psi); spsi=sin(psi);

rx = inf; ry=inf;
if abs(cpsi)>1e-12
    if cpsi>0, rx=xmax/cpsi; else, rx=xmin/cpsi; end
end
if abs(spsi)>1e-12
    if spsi>0, ry=ymax/spsi; else, ry=ymin/spsi; end
end

rmax = min([rx,ry]);
if rmax<0, rmax=0; end

end
