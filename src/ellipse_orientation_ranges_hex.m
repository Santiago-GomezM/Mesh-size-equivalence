function [out, int_def] = ellipse_orientation_ranges_hex(p, a, b, C)
% C es Nx2 con centros [cx_i, cy_i]

    N = size(C,1);

    % Prealocaciones
    out     = zeros(N,1);
    int_def = cell(N,1);

    % Normales del hexágono (fijas)
    phis = pi/6 + (0:5)*pi/3;   % 1x6
    ux   = cos(phis);           % 1x6
    uy   = sin(phis);           % 1x6

    % Parámetros de la elipse
    A2 = a*a;
    B2 = b*b;
    denom = A2 - B2;
    tol = 1e-12;

    % Precomputamos tabla de normales (para broadcasting)
    ux = reshape(ux,1,6);     % 1x6
    uy = reshape(uy,1,6);
    phis = reshape(phis,1,6);

    % Proyecciones d = ux*cx + uy*cy  → Nx6 totalmente vectorizado
    d = C(:,1).*ux + C(:,2).*uy;      % Nx6
    rhs = p - d;                      % Nx6

    % Casos en los que ya no cabe ninguna orientación
    invalid = any(rhs <= 0, 2);
    out(invalid) = 0;
    int_def(invalid) = {[]};

    % Si la elipse es (casi) circular
    if abs(denom) < 1e-15
        r = sqrt(A2);

        ok = ~invalid & all(rhs >= r - tol, 2);
        out(ok) = 1;
        int_def(ok) = { [0 180] };

        return
    end

    % K para *todos* los centros y caras → Nx6
    K = (2*rhs.^2 - (A2 + B2)) ./ denom;

    %-----------------------------------------------------------
    % LOOP sobre N, pero con toda la parte geométrica vectorizada
    %-----------------------------------------------------------
    for i = 1:N

        if invalid(i)
            continue
        end

        Ki   = K(i,:);       % 1x6
        rhsi = rhs(i,:);     % 1x6

        % Intervalos vigentes
        current = [0, pi];

        for k = 1:6

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

                th1 = phi - alpha2;
                th2 = phi - alpha1;

                a1 = mod(th1, pi);
                b1 = mod(th2, pi);

                if a1 <= b1
                    allowed = [a1, b1];
                else
                    allowed = [0, b1;
                               a1, pi];
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
            out(i) = sum(L)/pi;
            int_def{i} = current * 180/pi;
        end
    end
end

% function [out, int_def] = ellipse_orientation_ranges_hex_fast(p, a, b, c)
% % ellipse_orientation_ranges_hex_fast
% % Devuelve:
% %   out     = |Theta(c)| / pi    (fracción de orientaciones en [0,pi))
% %   int_def = Mx2 con intervalos en GRADOS [0,180]
% %
% % Parámetros:
% %   p  - apotema del hexágono regular centrado en el origen
% %   a,b - semiejes mayor y menor de la elipse
% %   c  - centro [cx,cy]
% 
%     cx = c(1); 
%     cy = c(2);
% 
%     % Normales del hexágono (ángulos de las caras)
%     phis = pi/6 + (0:5)*pi/3;  % 6x1
%     ux   = cos(phis);
%     uy   = sin(phis);
% 
%     % Proyección del centro sobre cada normal
%     d   = ux*cx + uy*cy;        % 6x1
%     rhs = p - d;                % 6x1
% 
%     % Tolerancia numérica
%     tol = 1e-12;
% 
%     % Si alguna cara ya está estrictamente "más allá" del apotema,
%     % no hay forma de que ninguna elipse quepa (independiente de theta).
%     % Esta condición es geométricamente correcta para elipse general.
%     if any(rhs <= 0)
%         out     = 0;
%         int_def = [];
%         return;
%     end
% 
%     A2 = a*a; 
%     B2 = b*b;
%     denom = A2 - B2;
% 
%     %---------------------------------------------------------
%     % CASO CIRCULAR / CASI CIRCULAR (a ≈ b)
%     %---------------------------------------------------------
%     if abs(denom) < 1e-15
%         % Elipse ~ círculo: soporte independiente de theta.
%         % El radio es r = a ~ b.
%         r = sqrt(A2);  % radio
% 
%         % Condición correcta de inclusión del círculo:
%         %   n_k · c + r <= p  <=>  rhs_k >= r  para todo k
%         %
%         % Para el límite r -> 0, cualquier punto interior debe dar
%         % rhs_k > 0 y, por tanto, todas las orientaciones son válidas.
%         %
%         % Usamos una pequeña tolerancia para robustez numérica.
%         if min(rhs) >= r - tol
%             out     = 1;             % todas las orientaciones en [0,pi)
%             int_def = [0 180];
%         else
%             out     = 0;
%             int_def = [];
%         end
%         return;
%     end
% 
%     %---------------------------------------------------------
%     % CASO ELÍPTICO GENERAL (a != b)
%     %---------------------------------------------------------
% 
%     % Parámetro K para cada cara
%     K = (2*rhs.^2 - (A2 + B2)) ./ denom;   % 6x1
% 
%     % Intervalos actuales de intersección en radianes
%     % Empezamos con todo [0,pi)
%     current = [0, pi];
% 
%     for k = 1:6
%         Kk  = K(k);
%         phi = phis(k);
% 
%         if Kk < -1
%             % Esta cara descarta todas las orientaciones
%             out     = 0;
%             int_def = [];
%             return;
% 
%         elseif Kk >= 1
%             % Esta cara no restringe nada: intersección no cambia
%             continue;
% 
%         else
%             % -1 <= Kk < 1 : hay un intervalo (o dos) de orientaciones válidas
% 
%             alpha0 = acos(Kk);        % en (0,pi)
%             alpha1 = 0.5*alpha0;      % límite inferior en α
%             alpha2 = pi - alpha1;     % límite superior en α
% 
%             % alpha = phi - theta ∈ [alpha1, alpha2]
%             % => theta ∈ [phi - alpha2, phi - alpha1]
%             th1 = phi - alpha2;
%             th2 = phi - alpha1;
% 
%             % Reducir a [0,pi) usando periodo pi
%             a1 = mod(th1, pi);
%             b1 = mod(th2, pi);
% 
%             if a1 <= b1
%                 % un único intervalo
%                 allowed = [a1, b1];
%             else
%                 % Envuelve: dos intervalos [0,b1] U [a1,pi]
%                 allowed = [0,  b1;
%                            a1, pi];
%             end
% 
%             % Intersección entre 'current' y 'allowed'
%             current = intersect_interval_sets(current, allowed);
% 
%             if isempty(current)
%                 out     = 0;
%                 int_def = [];
%                 return;
%             end
%         end
%     end
% 
%     % current contiene ahora Mx2 intervalos en radianes en [0,pi)
%     current = merge_intervals(current);
% 
%     lengths   = current(:,2) - current(:,1);
%     total_len = sum(lengths);
% 
%     out = total_len / pi;  % fracción de orientaciones en [0,pi)
% 
%     % Devolvemos intervalos en grados, como tu código original
%     int_def = current * 180/pi;
% end

%---------------------- AUXILIARES ---------------------------%

function S = intersect_interval_sets(A, B)
% A: Mx2, B: Nx2, todos en [0,pi)
% Devuelve la intersección como lista de intervalos (posiblemente vacía).

    if isempty(A) || isempty(B)
        S = [];
        return;
    end

    S = [];
    for i = 1:size(A,1)
        a1 = A(i,1); b1 = A(i,2);
        for j = 1:size(B,1)
            a2 = B(j,1); b2 = B(j,2);
            L = max(a1,a2);
            R = min(b1,b2);
            if R > L
                S = [S; L, R];
            end
        end
    end

    if isempty(S)
        return;
    end

    % Fusionar posibles solapamientos
    S = merge_intervals(S);
end

function I = merge_intervals(I)
% Fusiona intervalos solapados o contiguos.
    if isempty(I)
        return;
    end
    I = sortrows(I,1);
    res = I(1,:);
    for k = 2:size(I,1)
        a = I(k,1); b = I(k,2);
        if a <= res(end,2) + 1e-12  % solapan o se tocan (con tolerancia)
            res(end,2) = max(res(end,2), b);
        else
            res = [res; a, b];
        end
    end
    I = res;
end

