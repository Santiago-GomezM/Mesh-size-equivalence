function ratio = equivalent_opening_hex(Gamma, Delta, N, qstar, ...
        lmin, lmax, tol)
%--------------------------------------------------------------------------
% Computes the ratio l1/l2 for given shape ratios b/a and c/a
%
% The ellipsoid is defined as:
%   a = 1, b = Gamma, c = Delta
%
% INPUT:
%   Gamma      : b/a
%   Delta      : c/a
%   N          : number of attempts
%   qstar      : target cumulative probability
%   lmin,lmax  : search interval for l
%   tol        : tolerance (default 1e-6)
%
% OUTPUT:
%   ratio      : l1 / l2
%--------------------------------------------------------------------------

    if nargin < 9
        tol = 1e-6;
    end

    %---------------------------------------------
    % Fix scale (no loss of generality)
    %---------------------------------------------
    a = 1.0;
    b = Gamma;
    c = Delta;

    %---------------------------------------------
    % cumulative probabilities
    %---------------------------------------------
    q1 = @(l) 1 - (1 - prob_ellipsoid_hex_fast_opt(l,1,Gamma,Delta)).^N;
    %q1 = @(l) 1 - (1 - prob_ellipsoid_rect_fast_opt([-l/2 l/2 -1.3889*l/2 1.3889*l/2],1,Gamma,Delta)).^N;
    q2 = @(l) 1 - (1 - prob_ellipsoid_rect_fast_opt([-l/2 l/2 -l/2 l/2],1,Gamma,Delta)).^N;

    %---------------------------------------------
    % solve q1(l1) = qstar
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
    % solve q2(l2) = qstar
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
    % ratio
    %---------------------------------------------
    ratio = l1 / l2;
end