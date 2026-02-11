function ratio = equivalent_opening_rect(Gamma, Delta, N, qstar, ...
        lmin, lmax, tol, aspect)
%EQUIVALENT_OPENING_RECT Equivalent opening ratio between rectangular and square apertures.
%
%   ratio = EQUIVALENT_OPENING_RECT(Gamma, Delta, N, qstar, lmin, lmax, tol, aspect)
%   computes the ratio l1/l2 between a rectangular opening and a square
%   opening such that both yield the same cumulative passage probability
%   qstar after N independent presentation attempts.
%
%   The particle is modeled as a triaxial ellipsoid with normalized
%   semi-axes:
%       a = 1,  b = Gamma,  c = Delta.
%
%   The cumulative probability of passage after N attempts is:
%       q(l) = 1 - (1 - p(l))^N,
%   where p(l) is the single-attempt passage probability for a given
%   aperture size l.
%
%   The function solves:
%       q1(l1) = qstar   (rectangular aperture with aspect ratio 'aspect')
%       q2(l2) = qstar   (square aperture)
%   using a bisection scheme, and returns ratio = l1 / l2.
%
%   Inputs:
%     Gamma  - shape ratio b/a
%     Delta  - shape ratio c/a
%     N      - number of presentation attempts
%     qstar  - target cumulative passage probability
%     lmin   - lower bound of search interval for l
%     lmax   - upper bound of search interval for l
%     tol    - bisection tolerance (default: 1e-6)
%     aspect - rectangle aspect factor applied as [-l/2, l/2] × [-aspect*l/2, aspect*l/2]
%
%   Output:
%     ratio  - equivalent opening ratio l1 / l2
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

    if nargin < 9
        tol = 1e-6;
    end

    %---------------------------------------------
    % Fix scale (without loss of generality)
    %---------------------------------------------
    a = 1.0;
    b = Gamma;
    c = Delta;

    %---------------------------------------------
    % Cumulative passage probabilities
    %---------------------------------------------
    q1 = @(l) 1 - (1 - prob_ellipsoid_rect([-l/2 l/2 -aspect*l/2 aspect*l/2],1,Gamma,Delta)).^N;
    q2 = @(l) 1 - (1 - prob_ellipsoid_rect([-l/2 l/2 -l/2 l/2],1,Gamma,Delta)).^N;

    %---------------------------------------------
    % Solve q1(l1) = qstar (rectangular opening)
    %---------------------------------------------
    f1 = @(l) q1(l) - qstar;
    l1_lo = lmin; 
    l1_hi = lmax;

    if f1(l1_lo) > 0 || f1(l1_hi) < 0
        error('Opening 1: root not bracketed.');
    end

    while (l1_hi - l1_lo) > tol
        lm = 0.5*(l1_lo + l1_hi);
        if f1(lm) >= 0
            l1_hi = lm;
        else
            l1_lo = lm;
        end
    end
    l1 = 0.5*(l1_lo + l1_hi);

    %---------------------------------------------
    % Solve q2(l2) = qstar (square opening)
    %---------------------------------------------
    f2 = @(l) q2(l) - qstar;
    l2_lo = lmin; 
    l2_hi = lmax;

    if f2(l2_lo) > 0 || f2(l2_hi) < 0
        error('Opening 2: root not bracketed.');
    end

    while (l2_hi - l2_lo) > tol
        lm = 0.5*(l2_lo + l2_hi);
        if f2(lm) >= 0
            l2_hi = lm;
        else
            l2_lo = lm;
        end
    end
    l2 = 0.5*(l2_lo + l2_hi);

    %---------------------------------------------
    % Equivalent opening ratio
    %---------------------------------------------
    ratio = l1 / l2;
end