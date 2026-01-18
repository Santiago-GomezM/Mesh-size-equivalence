function [X, Y, percepass] = passing_area_Hexagon_intersection(p, a, b, nTheta, nPhi)
% passing_area_Hexagon_INTERSECTION
% Calcula la región donde la elipse E(a,b) cabe en el hexágono
% para **todas** las orientaciones theta ∈ [0, 2π].

if nargin < 4, nTheta = 500; end
if nargin < 5, nPhi = 500; end

alpha = pi/6 + (0:5)*pi/3;      % Normales del hexágono (6)
theta = linspace(0, 2*pi, nTheta);

% Precompute support function h_k(theta)
ang = alpha' - theta;           % 6 x nTheta
Delta = sqrt((a*cos(ang)).^2 + (b*sin(ang)).^2);   % h_k(theta)

phi = linspace(0, 2*pi, nPhi);
r = zeros(1, nPhi);

for j = 1:nPhi
    denom = abs(cos(phi(j) - alpha));   % 1 x 6
    denom(denom < 1e-12) = NaN;         % evitar 0
    denomMat = repmat(denom', 1, nTheta);

    aux = p - Delta;       % p - h_k(theta)
    aux(aux < 0) = 0;      % si ya no cabe, radio = 0

    limits = aux ./ denomMat;   % 6 x nTheta

    % *** INTERSECCIÓN ***
    % Mínimo en θ (peor orientación posible)
    minLimits = min(limits, [], 1);     % 1 x nTheta

    % Ahora el mínimo radial entre orientaciones
    r(j) = min(minLimits);              % <-- CAMBIO CLAVE
end

% A coordenadas cartesianas
X = r .* cos(phi);
Y = r .* sin(phi);

areaPermissive = polyarea(X, Y);
hexArea = 2*sqrt(3)*p^2;
percepass = areaPermissive / hexArea;

end