function varargout = wgs84Geodetic(varargin)
%WGS84GEODETIC  Convert between ECEF Cartesian and WGS84 geodetic coordinates.
%
%   ECEF Cartesian -> Geodetic:
%     [lat_deg, lon_deg, alt_km] = wgs84Geodetic(r_ecef_km)
%     where r_ecef_km is a 3-element vector (3x1 or 1x3).
%     If r_ecef_km is 3xN, outputs lat/lon/alt are 1xN row vectors.
%
%   Geodetic -> ECEF Cartesian:
%     r_ecef_km = wgs84Geodetic(lat_deg, lon_deg, alt_km, 'inverse')
%     Returns 3x1 vector.
%
%   WGS84 parameters:
%     a_E = 6378.137 km, f = 1/298.257223563
%     b_E = a_E*(1-f) = 6356.7523 km
%     e^2 = 1 - (b_E/a_E)^2 = 0.00669438 (first eccentricity squared)
%
%   ECEF -> Geodetic algorithm: Bowring iterative method (3-4 iterations).
%   Geodetic -> ECEF: closed-form using prime vertical radius N.

%% ── WGS84 constants ──────────────────────────────────────────────────────────
a_E = 6378.137;                % km, semi-major axis
f   = 1 / 298.257223563;
b_E = a_E * (1 - f);           % km, semi-minor axis
e2  = 1 - (b_E / a_E)^2;       % first eccentricity squared

%% ── Parse inputs to determine direction ─────────────────────────────────────
if nargin >= 2 && ischar(varargin{end}) && strcmpi(varargin{end}, 'inverse')
    %% ── Geodetic -> ECEF ─────────────────────────────────────────────────────
    if nargin < 4
        error('wgs84Geodetic: inverse mode requires (lat_deg, lon_deg, alt_km, ''inverse'')');
    end
    lat_deg = varargin{1};
    lon_deg = varargin{2};
    alt_km  = varargin{3};

    N_pv = a_E ./ sqrt(1 - e2 .* sind(lat_deg).^2);   % prime vertical radius

    r = [(N_pv + alt_km) .* cosd(lat_deg) .* cosd(lon_deg); ...
         (N_pv + alt_km) .* cosd(lat_deg) .* sind(lon_deg); ...
         (N_pv .* (1 - e2) + alt_km) .* sind(lat_deg)];

    varargout{1} = r;

else
    %% ── ECEF -> Geodetic ─────────────────────────────────────────────────────
    if nargin < 1
        error('wgs84Geodetic: requires at least r_ecef_km');
    end
    r_ecef = varargin{1};

    % Accept 3x1, 1x3, or 3xN
    if isvector(r_ecef) && numel(r_ecef) == 3
        r_ecef = r_ecef(:);   % ensure 3x1
        N_pts = 1;
    else
        % Expect 3xN
        if size(r_ecef, 1) ~= 3
            error('wgs84Geodetic: r_ecef_km must be a 3-element vector or 3xN matrix');
        end
        N_pts = size(r_ecef, 2);
    end

    lat_out = zeros(1, N_pts);
    lon_out = zeros(1, N_pts);
    alt_out = zeros(1, N_pts);

    for k = 1:N_pts
        r = r_ecef(:, k);

        % Distance from Z-axis
        p = sqrt(r(1)^2 + r(2)^2);

        % Longitude
        lon = atan2d(r(2), r(1));

        % Bowring iterative latitude
        % Initial estimate (geocentric latitude adjusted for ellipsoid)
        lat0 = atan2d(r(3), p * (1 - e2));

        for iter = 1:15
            N_pv   = a_E / sqrt(1 - e2 * sind(lat0)^2);
            lat_new = atan2d(r(3) + e2 * N_pv * sind(lat0), p);
            if abs(lat_new - lat0) < 1e-12
                break;
            end
            lat0 = lat_new;
        end

        N_pv = a_E / sqrt(1 - e2 * sind(lat0)^2);

        % Altitude: use p/cos(lat) formula, but near poles use z/sin(lat)
        if abs(lat0) > 89.0
            alt = r(3) / sind(lat0) - N_pv * (1 - e2);
        else
            alt = p / cosd(lat0) - N_pv;
        end

        lat_out(k) = lat0;
        lon_out(k) = lon;
        alt_out(k) = alt;
    end

    varargout{1} = lat_out;
    varargout{2} = lon_out;
    varargout{3} = alt_out;
end
end
