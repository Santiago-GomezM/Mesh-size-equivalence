function A = area_Dtheta_rect(rect, a, b, theta)
% AREA_DTHETA_RECT_CLIP
%   Calcula el área del conjunto D_theta = R ⊖ E_theta,
%   donde R es un rectángulo axis-aligned y E_theta es una elipse
%   de semiejes (a,b) rotada un ángulo theta.
%
% INPUT:
%   rect  : [xmin xmax ymin ymax]
%   a, b  : semiejes de la elipse (a >= b)
%   theta : orientación de la elipse (rad)
%
% OUTPUT:
%   A : área del polígono D_theta (0 si el conjunto es vacío)

    xmin = rect(1); xmax = rect(2);
    ymin = rect(3); ymax = rect(4);

    %--- 1) Normales del rectángulo (caras)
    normals = [ 1 0;   % x <= xmax
               -1 0;   % x >= xmin
                0 1;   % y <= ymax
                0 -1]; % y >= ymin

    d = [xmax; -xmin; ymax; -ymin];

    %--- 2) Matriz de la elipse rotada
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    M = R * diag([a*a, b*b]) * R.';

    %--- 3) Soporte h_E_theta(n_k)
    h = zeros(4,1);
    for k = 1:4
        u = normals(k,:)';
        h(k) = sqrt(u' * M * u);
    end

    %--- 4) Semiespacios n_k·C <= r_k(theta)
    r = d - h;

    %--- 5) Polígono inicial grande
    %     Tomamos un cuadrado suficientemente grande que seguro
    %     contiene D_theta si no es vacío.
    L = max([abs(xmin), abs(xmax), abs(ymin), abs(ymax)]) + max(a,b);
    P = [ -L -L;
           L -L;
           L  L;
          -L  L ];

    %--- 6) Clipping sucesivo por los 4 semiespacios
    for k = 1:4
        nk = normals(k,:);
        rk = r(k);
        P  = clip_polygon_halfspace(P, nk, rk);
        if isempty(P)
            A = 0;
            return;
        end
    end

    %--- 7) Área del polígono resultante
    A = polyarea(P(:,1), P(:,2));
end


% =============================================================
% Clipping de un polígono P por el semiespacio {x : n·x <= r}
% =============================================================
function Q = clip_polygon_halfspace(P, n, r)
    if isempty(P)
        Q = [];
        return;
    end

    Q = [];
    N = size(P,1);

    for i = 1:N
        A = P(i,:);
        B = P(mod(i,N)+1,:);

        fA = dot(n,A) - r;
        fB = dot(n,B) - r;

        insideA = (fA <= 0);
        insideB = (fB <= 0);

        if insideA && insideB
            Q = [Q; B];
        elseif insideA && ~insideB
            t = fA / (fA - fB);
            I = A + t*(B-A);
            Q = [Q; I];
        elseif ~insideA && insideB
            t = fA / (fA - fB);
            I = A + t*(B-A);
            Q = [Q; I; B];
        end
    end
end