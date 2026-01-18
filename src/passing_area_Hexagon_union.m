function [X, Y, percepass] = passing_area_Hexagon_union(p, a, b, nTheta, nPhi)
% passing_area_Hexagon_fast - Vectorized version

if nargin < 4, nTheta = 500; end
if nargin < 5, nPhi   = 500; end

alpha = pi/6 + (0:5)*pi/3;          % 1×6 normales del hexágono
theta = linspace(0, 2*pi, nTheta);  % 1×nTheta orientaciones de la elipse

%% ---- Precompute Delta (6 × nTheta) ----
ang   = alpha.' - theta;            % 6×nTheta
Delta = sqrt((a*cos(ang)).^2 + (b*sin(ang)).^2);

%% ---- Cálculo vectorizado para todas las phi ----
phi = linspace(0, 2*pi, nPhi);      % 1×nPhi

% OJO: aquí va ABS, como en tu código original
denom = abs(cos(phi - alpha.'));    % 6×nPhi
denom(denom < 1e-9) = NaN;          % evitar división por casi-cero

% aux: margen entre borde del hexágono y el de la elipse
aux = p - Delta;                    % 6×nTheta
aux(aux < 0) = 0;                   % no permitimos valores negativos

% Dar forma 3D para usar expansión implícita:
%   aux3 : 6×nTheta×1
%   den3 : 6×1×nPhi
aux3 = reshape(aux,   6, nTheta, 1);
den3 = reshape(denom, 6, 1,      nPhi);

% limits: 6 × nTheta × nPhi
limits = aux3 ./ den3;

% mínimo sobre los 6 lados del hexágono → nTheta×nPhi
minAcrossSides = squeeze(min(limits, [], 1));  % nTheta×nPhi

% máximo sobre las orientaciones de la elipse → 1×nPhi
r = max(minAcrossSides, [], 1);                % 1×nPhi

%% ---- Pasar a cartesiano ----
X = r .* cos(phi);
Y = r .* sin(phi);

areaPermissive = polyarea(X, Y);
hexArea        = 2*sqrt(3)*p^2;
percepass      = areaPermissive / hexArea;

end

% function [X, Y, percepass] = passing_area_Hexagon(p, a, b, nTheta, nPhi)
% % passing_area_Hexagon - Computes the permissive region for an ellipse inside a hexagon
% %
% % Inputs:
% %   p      - Apothem of the hexagon
% %   a, b   - Semi-major and semi-minor axes of the ellipse
% %   nTheta - Number of orientation samples (default: 500)
% %   nPhi   - Number of radial directions for boundary (default: 500)
% %
% % Outputs:
% %   X, Y        - Coordinates of the permissive region boundary
% %   percepass   - Ratio of permissive area to hexagon area
% 
% if nargin < 4, nTheta = 500; end
% if nargin < 5, nPhi = 500; end
% 
% alpha = pi/6 + (0:5)*pi/3;          % Hexagon normals
% theta = linspace(0, 2*pi, nTheta);  % Ellipse orientations
% 
% % Precompute Delta for all sides and orientations (vectorized)
% ang = alpha' - theta;  % 6 x nTheta matrix
% Delta = sqrt((a*cos(ang)).^2 + (b*sin(ang)).^2);
% 
% phi = linspace(0, 2*pi, nPhi);
% r = zeros(1, nPhi);
% 
% for j = 1:nPhi
%     denom = abs(cos(phi(j) - alpha)); % 1 x 6
%     denom(denom < 1e-9) = NaN;        % Avoid division by zero
%     % Expand denom to match Delta size (6 x nTheta)
%     denomMat = repmat(denom', 1, nTheta);
%     aux = (p - Delta);
%     aux(aux<0) = 0;
%     limits = aux ./ denomMat; % 6 x nTheta
%     minLimits = min(limits, [], 1);   % min across sides
%     r(j) = max(minLimits);            % best orientation for this direction
% end
% 
% % Convert to Cartesian
% X = r .* cos(phi);
% Y = r .* sin(phi);
% 
% areaPermissive = polyarea(X, Y);
% hexArea = 2*sqrt(3)*p^2;
% percepass = areaPermissive / hexArea;
% 
% end

% function [X, Y, percepass] = passing_area_Hexagon(p, a, b, nTheta, nPhi)
% % passing_area_Hexagon - Computes the permissive region for an ellipse inside a hexagon
% %
% % Syntax:
% %   [X, Y, percepass] = passing_area_Hexagon(p, a, b, nTheta, nPhi)
% %
% % Inputs:
% %   p      - Apothem of the hexagon
% %   a, b   - Semi-major and semi-minor axes of the ellipse
% %   nTheta - Number of orientation samples (default: 500)
% %   nPhi   - Number of radial directions for boundary (default: 500)
% %
% % Outputs:
% %   X, Y        - Coordinates of the permissive region boundary
% %   percepass   - Ratio of permissive area to hexagon area
% %
% % Example:
% %   [X, Y, ratio] = passing_area_Hexagon(5, 1, 0.8, 500, 500);
% %   figure; fill(X, Y, [0.7 0.9 1]); axis equal;
% 
% if nargin < 4, nTheta = 500; end
% if nargin < 5, nPhi = 500; end
% 
% % Hexagon normals
% alpha = pi/6 + (0:5)*pi/3; % angles of normals
% theta = linspace(0, 2*pi, nTheta); % orientations of ellipse
% 
% % Compute offsets for each orientation and each side
% Delta = zeros(6, nTheta);
% for i = 1:6
%     for k = 1:nTheta
%         ang = alpha(i) - theta(k);
%         Delta(i,k) = sqrt((a*cos(ang))^2 + (b*sin(ang))^2);
%     end
% end
% 
% % Compute boundary in polar coordinates
% phi = linspace(0, 2*pi, nPhi);
% r = zeros(size(phi));
% 
% for j = 1:nPhi
%     d_vals = zeros(1, nTheta);
%     for k = 1:nTheta
%         limits = Inf(1,6); % initialize with Inf
%         for i = 1:6
%             denom = abs(cos(phi(j) - alpha(i)));
%             if denom > 1e-9 % only consider positive denominators
%                 val = (p - Delta(i,k)) / denom;
%                 limits(i) = val;
%             end
%         end
%         % Take minimum of valid limits
%         minLimit = min(limits);
%         d_vals(k) = minLimit;
%     end
%     r(j) = max(d_vals); % best orientation for this direction
% end
% 
% % Convert to Cartesian
% X = r .* cos(phi);
% Y = r .* sin(phi);
% 
% % Compute area
% areaPermissive = polyarea(X, Y);
% hexArea = 2*sqrt(3)*p^2; % area of hexagon
% percepass = areaPermissive / hexArea;
% 
% end
% 
