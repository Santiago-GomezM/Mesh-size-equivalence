function [out, int_def] = ellipse_orientation_ranges_rect(rect,a,b,C)
% ELLIPSE_ORIENTATION_RANGES_RECT_FAST
%   Calcula los intervalos de orientación θ ∈ [0,π) para los que una
%   elipse de semiejes (a,b) centrada en C cabe dentro de un rectángulo
%   axis-aligned RECT = [xmin xmax ymin ymax].
%
% INPUT
%   a,b   : semiejes de la elipse (a >= b)
%   C     : Nx2, centros [cx_i, cy_i]
%   rect  : [xmin xmax ymin ymax]
%
% OUTPUT
%   out     : Nx1, fracción de orientaciones viables |Θ_i| / π
%   int_def : cell{N,1}, intervalos viables en grados

    N = size(C,1);

    %--------------------------------------------------
    % Prealocación
    %--------------------------------------------------
    out     = zeros(N,1);
    int_def = cell(N,1);

    %--------------------------------------------------
    % Geometría fija del rectángulo
    %--------------------------------------------------
    normals = [ 1  0;    % x <= xmax
               -1  0;    % x >= xmin
                0  1;    % y <= ymax
                0 -1 ];  % y >= ymin

    d = [rect(2); -rect(1); rect(4); -rect(3)];   % términos independientes

    phis = [0; pi; pi/2; -pi/2];                  % ángulos de las normales

    %--------------------------------------------------
    % Parámetros de la elipse
    %--------------------------------------------------
    A2 = a*a;
    B2 = b*b;
    denom = A2 - B2;
    tol = 1e-12;

    %--------------------------------------------------
    % Parte geométrica vectorizada
    %--------------------------------------------------
    % r(i,k) = margen disponible del centro i frente a la cara k
    r = d.' - C*normals.';        % Nx4

    % Centros inválidos (ni la elipse mínima cabe)
    invalid = any(r <= 0, 2);
    out(invalid) = 0;
    int_def(invalid) = {[]};

    %--------------------------------------------------
    % Caso casi circular
    %--------------------------------------------------
    if abs(denom) < 1e-15
        h = sqrt(A2);
        ok = ~invalid & all(r >= h - tol, 2);
        out(ok) = 1;
        int_def(ok) = { [0 180] };
        return
    end

    %--------------------------------------------------
    % Coeficientes K para todos los centros y caras
    %--------------------------------------------------
    K = (2*r.^2 - (A2 + B2)) ./ denom;   % Nx4

    %--------------------------------------------------
    % Loop lógico (intervalos) sobre centros
    %--------------------------------------------------
    for i = 1:N

        if invalid(i)
            continue
        end

        Ki = K(i,:);    % 1x4

        % Intervalos iniciales: todo [0,π)
        current = [0 pi];

        for k = 1:4

            Kk  = Ki(k);
            phi = phis(k);

            if Kk < -1
                current = [];
                break

            elseif Kk >= 1
                continue

            else
                alpha0 = acos(Kk);
                alpha1 = 0.5*alpha0;
                alpha2 = pi - alpha1;

                % theta ∈ [phi - alpha2, phi - alpha1]
                th1 = phi - alpha2;
                th2 = phi - alpha1;

                a1 = mod(th1, pi);
                b1 = mod(th2, pi);

                if a1 <= b1
                    allowed = [a1 b1];
                else
                    allowed = [0 b1;
                               a1 pi];
                end

                current = intersect_interval_sets(current, allowed);
                if isempty(current)
                    break
                end
            end
        end

        if isempty(current)
            out(i) = 0;
            int_def{i} = [];
        else
            current = merge_intervals(current);
            L = current(:,2) - current(:,1);
            out(i) = sum(L) / pi;
            int_def{i} = current * 180/pi;
        end
    end
end

%====================================================================
% AUXILIARES DE INTERVALOS EN [0,π)
%====================================================================
function S = intersect_interval_sets(A,B)
    if isempty(A) || isempty(B)
        S = [];
        return
    end
    S = [];
    for i = 1:size(A,1)
        for j = 1:size(B,1)
            L = max(A(i,1), B(j,1));
            R = min(A(i,2), B(j,2));
            if R > L
                S = [S; L R]; %#ok<AGROW>
            end
        end
    end
    S = merge_intervals(S);
end

function I = merge_intervals(I)
    if isempty(I), return; end
    I = sortrows(I,1);
    res = I(1,:);
    for k = 2:size(I,1)
        if I(k,1) <= res(end,2) + 1e-12
            res(end,2) = max(res(end,2), I(k,2));
        else
            res = [res; I(k,:)]; %#ok<AGROW>
        end
    end
    I = res;
end