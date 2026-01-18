function [l,N,it] = equivalent_square(a,b,c,l0,N0,eps)

if nargin < 6, eps = 1e-3; end

tol = 1e-4;
maxIter = 20;
delta = 1e-4;

l = l0;
N = N0;

for it = 1:maxIter

    % Probabilidades básicas
    p1 = prob_ellipsoid_rect_fast_opt([-l/2 l/2 -l/2 l/2], a,b,c);
    p2 = prob_ellipsoid_rect_fast_opt([-l l -l l], a,b,c);

    q1 = 1 - p1;
    q2 = 1 - p2;

    % Funciones
    F1 = 1 - q1^N - 0.9;
    F2 = 1 - q2^N - (1 - eps);

    F = [F1; F2];

    if norm(F,inf) < tol
        return
    end

    % Derivadas p'(l)
    p1p = prob_ellipsoid_rect_fast_opt([-(l+delta)/2 (l+delta)/2 ...
                                        -(l+delta)/2 (l+delta)/2], a,b,c);
    p1m = prob_ellipsoid_rect_fast_opt([-(l-delta)/2 (l-delta)/2 ...
                                        -(l-delta)/2 (l-delta)/2], a,b,c);
    dp1 = (p1p - p1m)/(2*delta);

    p2p = prob_ellipsoid_rect_fast_opt([-(2*l+delta) (2*l+delta) ...
                                        -(2*l+delta) (2*l+delta)], a,b,c);
    p2m = prob_ellipsoid_rect_fast_opt([-(2*l-delta) (2*l-delta) ...
                                        -(2*l-delta) (2*l-delta)], a,b,c);
    dp2 = (p2p - p2m)/(2*delta);

    % Jacobiano
    J = [ N*q1^(N-1)*dp1,   -q1^N*log(q1)
          2*N*q2^(N-1)*dp2, -q2^N*log(q2) ];

    % Paso de Newton
    deltaX = -J\F;
    l = l + deltaX(1);
    N = max(1, N + deltaX(2));   % N > 0
end

error('Newton no convergió');
end