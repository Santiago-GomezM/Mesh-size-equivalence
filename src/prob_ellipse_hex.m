function P = prob_ellipse_hex(p, a, b)

A_H = 2*sqrt(3)*p^2;   % área del hexágono regular de apotema p
g = @(theta) 1/pi;     % orientación uniforme
integrand = @(theta) g(theta) .* arrayfun(@(t) area_Dtheta_hex(p,a,b,t), theta);
P = (1/A_H) * integral(integrand, 0, pi, 'AbsTol',1e-8,'RelTol',1e-6);

end