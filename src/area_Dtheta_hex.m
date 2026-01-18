function A = area_Dtheta_hex(p, a, b, theta)
% AREA_DTHETA_HEX_CLIP
%   Devuelve el área del polígono D_theta = H ⊖ E_theta,
%   donde H es un hexágono regular centrado en el origen de apotema p,
%   y E_theta es una elipse de semiejes (a,b) rotada un ángulo theta.
%
% INPUT:
%   p     - apotema del hexágono
%   a, b  - semiejes de la elipse (a = semieje mayor, b = menor)
%   theta - orientación de la elipse (radianes)
%
% OUTPUT:
%   A - área de D_theta (0 si el conjunto es vacío)

    %--- 1) Normales del hexágono (6 caras) ---
    phis    = pi/6 + (0:5)'*pi/3;
    normals = [cos(phis), sin(phis)];   % 6 x 2

    %--- 2) Matriz de la elipse rotada ---
    R = [cos(theta) -sin(theta); 
         sin(theta)  cos(theta)];
    M = R * diag([a^2, b^2]) * R.';     % 2x2

    %--- 3) Soportes h_{E_theta}(n_k) ---
    h = zeros(6,1);
    for k = 1:6
        u = normals(k,:).';            % columna 2x1
        h(k) = sqrt(u' * M * u);
    end

    %--- 4) Semiespacios n_k·C <= r_k(theta) ---
    r = p - h;

    % Nota IMPORTANTE:
    % NO hacemos ningún "if any(r<=0) => vacío", porque eso NO es correcto
    % en general: puede haber región no vacía con algunos r_k < 0.
    % Dejamos que el algoritmo de recorte decida.

    %--- 5) Polígono inicial grande que contiene toda la posible D_theta ---
    % Radio circunscrito del hexágono:
    R_hex = 2*p/sqrt(3);
    L = R_hex + max(a,b);   % margen de seguridad
    P = [ -L -L;
           L -L;
           L  L;
          -L  L ];          % cuadrado grande

    %--- 6) Recorte sucesivo por cada semiespacio n_k·x <= r_k ---
    for k = 1:6
        nk = normals(k,:);
        rk = r(k);
        P = clip_polygon_halfspace(P, nk, rk);
        if isempty(P)
            A = 0;
            return;
        end
    end

    %--- 7) Área del polígono resultante ---
    A = polyarea(P(:,1), P(:,2));
end

%==================================================================
% Recorte de un polígono P por el semiespacio { x : n·x <= r }
%==================================================================

function Q = clip_polygon_halfspace(P, n, r)
    % P: Nx2, vértices en orden CCW o CW
    % n: 1x2, vector normal
    % r: escalar
    %
    % Devuelve Q: Mx2, polígono recortado (puede ser vacío).

    if isempty(P)
        Q = [];
        return;
    end

    Q = [];
    N = size(P,1);

    for i = 1:N
        A = P(i,:);
        B = P(mod(i,N)+1,:);   % siguiente vértice (cierre cíclico)

        fA = dot(n,A) - r;
        fB = dot(n,B) - r;

        % Clasificación de los puntos respecto al semiespacio
        insideA = (fA <= 0);
        insideB = (fB <= 0);

        if insideA && insideB
            % Ambos dentro: mantenemos B
            Q = [Q; B];

        elseif insideA && ~insideB
            % A dentro, B fuera: añadimos sólo el punto de intersección
            t = fA / (fA - fB);
            I = A + t*(B - A);
            Q = [Q; I];

        elseif ~insideA && insideB
            % A fuera, B dentro: primero el punto de corte, luego B
            t = fA / (fA - fB);
            I = A + t*(B - A);
            Q = [Q; I; B];

        else
            % A y B fuera: no añadimos nada
        end
    end

    % Q queda vacío si el polígono ha sido completamente recortado
end