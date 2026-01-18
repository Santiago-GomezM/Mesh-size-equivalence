function [X, Y, rI] = passing_area_Rectangle_intersection(a,b,rect,nTheta,nPsi)
% PASSING_AREA_RECT_INTERSECTION
%   Área de centros que permiten paso de la elipse para TODAS las
%   orientaciones (intersección sobre theta).
%
%   a,b  : semiejes de la elipse
%   rect : [xmin xmax ymin ymax]
%   nTheta : nº de muestras en orientación (default 400)
%   nPsi   : nº de direcciones angulares en polar (default 720)
%
%   Devuelve:
%     X,Y : coordenadas cartesianas de la frontera en polar
%     rI  : radios r_cap(psi) en cada dirección psi

    if nargin < 4, nTheta = 400; end
    if nargin < 5, nPsi   = 720; end

    % Muestra de ángulos polares
    psis   = linspace(-pi, pi, nPsi);
    thetas = linspace(0,  pi, nTheta);

    % Normales del rectángulo y términos independientes
    xmin=rect(1); xmax=rect(2); ymin=rect(3); ymax=rect(4);
    normals = [1 0; -1 0; 0 1; 0 -1];
    d       = [xmax; -xmin; ymax; -ymin];

    A2 = a*a; B2 = b*b;

    rI = zeros(size(psis));

    for j = 1:nPsi
        psi = psis(j);
        % Radio máximo dentro del rectángulo en dirección psi
        rmax = boundary_rect_polar(rect, psi);

        % Si la dirección no entra en el rectángulo, radio cero
        if rmax <= 0
            rI(j) = 0;
            continue;
        end

        % Queremos r_cap(psi) = inf_theta r_theta(psi)
        rmin = inf;

        for t = 1:nTheta
            th = thetas(t);

            % Soportes de la elipse en las 4 normales
            h = zeros(4,1);
            for k = 1:4
                phi_k = atan2(normals(k,2), normals(k,1));
                h(k) = sqrt(A2*cos(phi_k - th)^2 + B2*sin(phi_k - th)^2);
            end

            % Para cada cara, el límite radial (si la normal "corta" la dirección)
            r_k = inf(4,1);
            for k = 1:4
                phi_k = atan2(normals(k,2), normals(k,1));
                cp = cos(psi - phi_k);
                if cp > 1e-12
                    r_k(k) = (d(k) - h(k)) / cp;
                end
            end

            % Si alguna cara no admite radio positivo, no hay r_theta > 0
            if any(r_k <= 0)
                r_theta = 0;
            else
                r_theta = min(r_k);
            end

            rmin = min(rmin, r_theta);
        end

        % Recortar con la frontera del rectángulo
        rI(j) = max(0, min(rmin, rmax));
    end

    X = rI .* cos(psis);
    Y = rI .* sin(psis);
end

%====================================================================
% Radio hasta la frontera del rectángulo en la dirección psi
%====================================================================
function rmax = boundary_rect_polar(rect, psi)
    xmin=rect(1); xmax=rect(2); ymin=rect(3); ymax=rect(4);
    cpsi = cos(psi); spsi = sin(psi);

    candidates = [];

    % Intersección con x = xmin y x = xmax
    if abs(cpsi) > 1e-12
        r1 = xmin / cpsi;
        r2 = xmax / cpsi;
        candidates = [candidates, r1, r2];
    end

    % Intersección con y = ymin y y = ymax
    if abs(spsi) > 1e-12
        r3 = ymin / spsi;
        r4 = ymax / spsi;
        candidates = [candidates, r3, r4];
    end

    % Nos quedamos con los radios positivos cuyo punto está realmente en el rectángulo
    r_valid = [];
    for r = candidates
        if r > 0
            x = r*cpsi; y = r*spsi;
            if x >= xmin-1e-9 && x <= xmax+1e-9 && ...
               y >= ymin-1e-9 && y <= ymax+1e-9
                r_valid(end+1) = r; %#ok<AGROW>
            end
        end
    end

    if isempty(r_valid)
        rmax = 0;
    else
        rmax = min(r_valid);
    end
end