function A = area_Dtheta_circ(R, a)
% AREA_DTHETA_CIRC
%   Área del conjunto D_theta = B_R ⊖ E_theta
%   donde B_R es un círculo de radio R y
%   E_theta es una elipse de semiejes (a,b).
%
%   NOTA: el resultado es independiente de theta y b.

    if R <= a
        A = 0;
    else
        A = pi * (R - a)^2;
    end
end