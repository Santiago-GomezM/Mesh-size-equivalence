function [l,N,it] = equivalent_square(a,b,c,l0,N0,eps)
%EQUIVALENT_SQUARE estimation of an equivalent square opening and presentation attempts
%    for estimating the maximum size of a fragment size distribution.
%
%   [l, N, it] = EQUIVALENT_SQUARE(a, b, c, l0, N0, eps) estimates an
%   equivalent square opening size l and an effective number of presentation
%   attempts N by solving a two-equation nonlinear system using Newton's method.
%
%   The system is built from cumulative non-passage probabilities based on the
%   single-attempt passage probability p(l) for an ellipsoidal particle with
%   semi-axes (a,b,c) through a square aperture. The probability model is
%   evaluated via PROB_ELLIPSOID_RECT_FAST_OPT.
%
%   Specifically, with q(l) = 1 - p(l), the method enforces two target
%   cumulative passage levels (as implemented below) and iteratively updates
%   (l, N) until convergence in the infinity norm.
%
%   Inputs:
%     a, b, c - ellipsoid semi-axes
%     l0      - initial guess for the equivalent square size l
%     N0      - initial guess for the effective number of attempts N
%     eps     - tolerance parameter for the second target level (default: 1e-3)
%
%   Outputs:
%     l       - estimated equivalent square opening size
%     N       - estimated effective number of presentation attempts (N > 0)
%     it      - number of Newton iterations performed
%
%   Notes:
%   - Derivatives dp/dl are approximated by centered finite differences.
%   - The iteration stops when ||F||_inf < tol or after maxIter steps.
%
%   Part of the code accompanying the manuscript:
%   "A geometric–probabilistic framework for size equivalence in mineral screening".
%
%   Copyright (c) 2026 Santiago Gómez et al.
%   Licensed under the MIT License (see LICENSE file).
%   If you use this code in academic work, please cite the associated publication.

if nargin < 6, eps = 1e-3; end

tol = 1e-4;
maxIter = 20;
delta = 1e-4;

l = l0;
N = N0;

for it = 1:maxIter

    % Single-attempt passage probabilities
    p1 = prob_ellipsoid_rect([-l/2 l/2 -l/2 l/2], a,b,c);
    p2 = prob_ellipsoid_rect([-l l -l l], a,b,c);

    q1 = 1 - p1;
    q2 = 1 - p2;

    % Nonlinear system (target cumulative passage levels)
    F1 = 1 - q1^N - 0.9;
    F2 = 1 - q2^N - (1 - eps);

    F = [F1; F2];

    if norm(F,inf) < tol
        return
    end

    % Finite-difference derivatives dp/dl
    p1p = prob_ellipsoid_rect([-(l+delta)/2 (l+delta)/2 ...
                                        -(l+delta)/2 (l+delta)/2], a,b,c);
    p1m = prob_ellipsoid_rect([-(l-delta)/2 (l-delta)/2 ...
                                        -(l-delta)/2 (l-delta)/2], a,b,c);
    dp1 = (p1p - p1m)/(2*delta);

    p2p = prob_ellipsoid_rect([-(2*l+delta) (2*l+delta) ...
                                        -(2*l+delta) (2*l+delta)], a,b,c);
    p2m = prob_ellipsoid_rect([-(2*l-delta) (2*l-delta) ...
                                        -(2*l-delta) (2*l-delta)], a,b,c);
    dp2 = (p2p - p2m)/(2*delta);

    % Jacobian matrix
    J = [ N*q1^(N-1)*dp1,   -q1^N*log(q1)
          2*N*q2^(N-1)*dp2, -q2^N*log(q2) ];

    % Newton step
    deltaX = -J\F;
    l = l + deltaX(1);
    N = max(1, N + deltaX(2));   % enforce N > 0
end

error('Newton did not converge.');

end
