function P = prob_ellipsoid_rect(rect,a,b,c)
%PROB_ELLIPSOID_RECT Single-attempt passage probability for a triaxial ellipsoid through a rectangular aperture.
%
%   P = PROB_ELLIPSOID_RECT(rect, a, b, c) computes the single-attempt passage
%   probability of a triaxial ellipsoid with semi-axes (a,b,c) through an
%   axis-aligned rectangular opening defined by rect = [xmin xmax ymin ymax].
%
%   The probability is obtained by averaging the admissible centroid area
%   over random particle orientations in SO(3). Orientations are parameterized
%   using ZYZ Euler angles and the Haar measure:
%
%       (1/pi^2) * sin(beta) d(alpha) d(beta) d(gamma),
%
%   with alpha,gamma ∈ [0, pi) and beta ∈ [0, pi/2] (using symmetry).
%
%   Implementation highlights:
%   - Uses explicit quadrature on SO(3) (Gauss–Legendre in beta, uniform trapezoid in
%     alpha and gamma).
%   - Computes the exact projected ellipse from the rotated ellipsoid quadratic form
%     (equivalent to Q = R*Dinv*R').
%   - Evaluates the admissible area via polygon clipping for a rectangle
%     (see AREA_DTHETA_RECT below).
%   - Self-contained: helper routines are included at the end of this file.
%
%   Inputs:
%     rect - rectangle bounds [xmin xmax ymin ymax]
%     a    - ellipsoid semi-axis (principal axis 1)
%     b    - ellipsoid semi-axis (principal axis 2)
%     c    - ellipsoid semi-axis (principal axis 3)
%
%   Output:
%     P    - single-attempt passage probability
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

% -------------------------------------------------
% Rectangle area
% -------------------------------------------------
A_rect = (rect(2)-rect(1))*(rect(4)-rect(3));

% -------------------------------------------------
% Quadrature on SO(3)
% -------------------------------------------------
Na = 20;   % alpha
Ng = 20;   % gamma
Nb = 14;   % beta (Gauss–Legendre)

alpha = (0:Na-1)*(pi/Na);  da = pi/Na;
gamma = (0:Ng-1)*(pi/Ng);  dg = pi/Ng;

ca = cos(alpha); sa = sin(alpha);
cg = cos(gamma); sg = sin(gamma);

[beta, wb] = gauss_legendre(Nb,0,pi/2);
cb = cos(beta); sb = sin(beta);

% Inverse squared semi-axes
d1 = 1/a^2; 
d2 = 1/b^2; 
d3 = 1/c^2;

Pacc = 0;

% =================================================
% Main loop
% =================================================
for ib = 1:Nb
    wbeta = wb(ib)*sb(ib);   % Haar weight: sin(beta) d beta
    cb_i = cb(ib); 
    sb_i = sb(ib);

    for ia = 1:Na
        ca_i = ca(ia); 
        sa_i = sa(ia);

        ca_cb = ca_i*cb_i;
        sa_cb = sa_i*cb_i;
        ca_sb = ca_i*sb_i;
        sa_sb = sa_i*sb_i;

        for ig = 1:Ng
            cg_i = cg(ig); 
            sg_i = sg(ig);

            % -----------------------------------------
            % Columns of R = Rz(alpha)*Ry(beta)*Rz(gamma) (ZYZ)
            % -----------------------------------------
            % r1
            r11 = ca_cb*cg_i - sa_i*sg_i;
            r21 = sa_cb*cg_i + ca_i*sg_i;
            r31 = -sb_i*cg_i;

            % r2
            r12 = -ca_cb*sg_i - sa_i*cg_i;
            r22 = -sa_cb*sg_i + ca_i*cg_i;
            r32 =  sb_i*sg_i;

            % r3
            r13 =  ca_sb;
            r23 =  sa_sb;
            r33 =  cb_i;

            % -----------------------------------------
            % Q2, q, Q33 from Q = R*Dinv*R'
            % -----------------------------------------
            Q11 = d1*r11*r11 + d2*r12*r12 + d3*r13*r13;
            Q22 = d1*r21*r21 + d2*r22*r22 + d3*r23*r23;
            Q12 = d1*r11*r21 + d2*r12*r22 + d3*r13*r23;

            q1  = d1*r11*r31 + d2*r12*r32 + d3*r13*r33;
            q2  = d1*r21*r31 + d2*r22*r32 + d3*r23*r33;

            Q33 = d1*r31*r31 + d2*r32*r32 + d3*r33*r33;

            invQ33 = 1/Q33;

            m11 = Q11 - (q1*q1)*invQ33;
            m22 = Q22 - (q2*q2)*invQ33;
            m12 = Q12 - (q1*q2)*invQ33;

            % -----------------------------------------
            % Projected ellipse parameters
            % -----------------------------------------
            tr   = m11 + m22;
            disc = hypot(m11-m22,2*m12);
            lam1 = 0.5*(tr-disc);
            lam2 = 0.5*(tr+disc);

            lam1 = max(lam1,realmin('double'));
            lam2 = max(lam2,realmin('double'));

            a_proj = 1/sqrt(lam1);
            b_proj = 1/sqrt(lam2);
            theta  = 0.5*atan2(2*m12,(m11-m22));

            % -----------------------------------------
            % Admissible area via rectangular clipping
            % -----------------------------------------
            A = area_Dtheta_rect(rect,a_proj,b_proj,theta);

            Pacc = Pacc + wbeta*A;
        end
    end
end

% -------------------------------------------------
% Final normalization
% -------------------------------------------------
P = (1/A_rect)*(1/pi^2)*(da*dg)*Pacc;

% =================================================
%                      SUBFUNCTIONS (SELF-CONTAINED)
% =================================================

function A = area_Dtheta_rect(rect,aE,bE,theta)
%AREA_DTHETA_RECT Admissible centroid area for a rectangular aperture (clipping implementation).

    xmin = rect(1); xmax = rect(2);
    ymin = rect(3); ymax = rect(4);

    normals = [ 1 0;
               -1 0;
                0 1;
                0 -1];

    d = [xmax; -xmin; ymax; -ymin];

    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
    M = R*diag([aE^2,bE^2])*R.';

    h = zeros(4,1);
    for k = 1:4
        u = normals(k,:)';
        h(k) = sqrt(u'*M*u);
    end

    r = d - h;

    L = max([abs(xmin),abs(xmax),abs(ymin),abs(ymax)]) + max(aE,bE);
    P0 = [-L -L; L -L; L L; -L L];

    for k = 1:4
        P0 = clip_polygon_halfspace(P0,normals(k,:),r(k));
        if isempty(P0)
            A = 0;
            return;
        end
    end

    A = polyarea(P0(:,1),P0(:,2));
end

function Q = clip_polygon_halfspace(P,n,r)
%CLIP_POLYGON_HALFSPACE Clip polygon P by the half-space { x : n·x <= r }.

    Q = [];
    Np = size(P,1);
    for i = 1:Np
        A = P(i,:);
        B = P(mod(i,Np)+1,:);
        fA = dot(n,A)-r;
        fB = dot(n,B)-r;
        insideA = (fA<=0);
        insideB = (fB<=0);
        if insideA && insideB
            Q = [Q; B];
        elseif insideA && ~insideB
            t = fA/(fA-fB);
            Q = [Q; A+t*(B-A)];
        elseif ~insideA && insideB
            t = fA/(fA-fB);
            Q = [Q; A+t*(B-A); B];
        end
    end
end

function [x,w] = gauss_legendre(n,aGL,bGL)
%GAUSS_LEGENDRE Gauss–Legendre nodes and weights on [aGL, bGL] (Golub–Welsch).

    i = (1:n-1)';
    betaGW = i./sqrt(4*i.^2-1);
    T = diag(betaGW,1)+diag(betaGW,-1);
    [V,D] = eig(T);
    x0 = diag(D);
    [x0,idx] = sort(x0);
    V = V(:,idx);
    w0 = 2*(V(1,:)').^2;
    x = (bGL-aGL)/2*x0 + (aGL+bGL)/2;
    w = (bGL-aGL)/2*w0;
end

end