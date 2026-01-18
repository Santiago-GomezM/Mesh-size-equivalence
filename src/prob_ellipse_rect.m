function P = prob_ellipse_rect(a,b,rect)

g = @(theta) 1/pi;            % orientación uniforme
A_rect = (rect(2)-rect(1)) * (rect(4)-rect(3));
integrand = @(th) g(th) .* arrayfun(@(t) area_Dtheta_rect(rect,a,b,t), th);
P = (1/A_rect) * integral(integrand, 0, pi, 'AbsTol',1e-8,'RelTol',1e-6);

end