function [X,Y,rU] = passing_area_Rectangle_union(a,b,rect,nTheta,nPsi)

if nargin<4, nTheta=400; end
if nargin<5, nPsi=720; end

phis = linspace(-pi,pi,nPsi);
thetas = linspace(0,pi,nTheta);

xmin=rect(1); xmax=rect(2); ymin=rect(3); ymax=rect(4);

% Normales del rectángulo
normals = [1 0; -1 0; 0 1; 0 -1];
d = [xmax; -xmin; ymax; -ymin];

A2=a*a; B2=b*b;

rU = zeros(size(phis));

for j=1:nPsi
    psi=phis(j);
    rmax=boundary_rect_polar(rect,psi);
    best_r=0;

    for t=1:nTheta
        th=thetas(t);

        h = zeros(4,1);
        for k=1:4
            phi_k = atan2(normals(k,2), normals(k,1));
            h(k)=sqrt(A2*cos(phi_k-th)^2 + B2*sin(phi_k-th)^2);
        end

        r_k = (d - h) ./ max(cos(psi - atan2(normals(:,2), normals(:,1))), 1e-12);

        r_k(r_k<=0) = 0;
        r_theta = min(r_k);

        best_r = max(best_r, r_theta);
    end

    rU(j) = min(best_r, rmax);
end

X = rU .* cos(phis);
Y = rU .* sin(phis);
end

% frontera polar del rectángulo
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


% function [X, Y, percepass] = passing_area_Rectangle(L, l, a, b, nTheta, nPhi)
%     % calcularRegionPermisiva: Calcula la región permisiva para una elipse dentro de un rectángulo
%     %
%     % Parámetros:
%     %   L      - Longitud del rectángulo
%     %   l      - Altura del rectángulo
%     %   a, b   - Semiejes de la elipse
%     %   nTheta - Número de orientaciones
%     %   nPhi   - Número de direcciones radiales
%     %
%     % Salidas:
%     %   areaPermisiva - Área aproximada de la región permisiva
%     %   percepass     - Porcentaje respecto al área del rectángulo
% 
%     % Orientaciones posibles
%     theta = linspace(0, pi/2, nTheta);
%     w = a*cos(theta) + b*sin(theta);
%     h = a*sin(theta) + b*cos(theta);
% 
%     % Ángulos para la frontera (direcciones radiales)
%     phi = linspace(0, 2*pi, nPhi);
%     r = zeros(size(phi));
% 
%     % Calcular el radio máximo para cada dirección phi
%     for j = 1:nPhi
%         d_vals = zeros(1, nTheta);
%         for i = 1:nTheta
%             % Evitar división por cero
%             if abs(cos(phi(j))) < 1e-9
%                 dx = Inf;
%             else
%                 dx = (L/2 - w(i)) / abs(cos(phi(j)));
%                 if dx < 0
%                     dx = 0;
%                 end
%             end
%             if abs(sin(phi(j))) < 1e-9
%                 dy = Inf;
%             else
%                 dy = (l/2 - h(i)) / abs(sin(phi(j)));
%                 if dy < 0
%                     dy = 0;
%                 end
%             end
%             d_vals(i) = min(dx, dy);
%         end
%         r(j) = max(d_vals);
%     end
% 
%     % Convertir a coordenadas cartesianas
%     X = r .* cos(phi);
%     Y = r .* sin(phi);
% 
%     % Calcular área aproximada
%     areaPermisiva = polyarea(X, Y);
%     percepass = areaPermisiva / (L * l);
% 
% end