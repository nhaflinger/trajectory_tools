function [az, el, range_km] = topocentricAzEl(obs_lat, obs_lon, obs_alt_km, r_sat_eci, jd)
%TOPOCENTRICAZEL  Compute azimuth and elevation from a ground observer to a satellite.
%
%   [az, el, range_km] = topocentricAzEl(obs_lat, obs_lon, obs_alt_km, r_sat_eci, jd)
%
%   For multiple time steps, r_sat_eci is 3xN and jd is 1xN (or scalar).
%
%   Inputs:
%     obs_lat     - observer geodetic latitude (deg, +North)
%     obs_lon     - observer longitude (deg, +East)
%     obs_alt_km  - observer altitude above mean sea level (km)
%     r_sat_eci   - satellite ECI position (km), 3xN or 3x1
%     jd          - Julian Date, 1xN or scalar
%
%   Outputs:
%     az          - azimuth (deg, 0=North, clockwise), 1xN
%     el          - elevation (deg, 0=horizon, +up), 1xN
%     range_km    - slant range to satellite (km), 1xN
%
%   Algorithm:
%     1. Compute observer ECEF (spherical approximation: R_E + alt)
%     2. Rotate satellite ECI → ECEF using GMST
%     3. Form range vector in ECEF, project onto ENU frame
%     4. Convert ENU to az/el

R_E = 6378.1363;   % km

%% ── Handle vectorized inputs ─────────────────────────────────────────────────
N = size(r_sat_eci, 2);   % number of time steps

% Expand scalar JD to match N
if isscalar(jd)
    jd = jd * ones(1, N);
end
jd = jd(:)';   % ensure row vector 1xN

%% ── Observer ECEF position (spherical Earth approximation) ───────────────────
r_obs_mag = R_E + obs_alt_km;
r_obs_ecef = r_obs_mag * [cosd(obs_lat)*cosd(obs_lon); ...
                           cosd(obs_lat)*sind(obs_lon); ...
                           sind(obs_lat)];   % 3x1, constant

%% ── ENU unit vectors at observer location (ECEF, constant) ───────────────────
e_E_ecef = [-sind(obs_lon);           cosd(obs_lon);          0           ];
e_N_ecef = [-sind(obs_lat)*cosd(obs_lon); -sind(obs_lat)*sind(obs_lon); cosd(obs_lat)];
e_U_ecef = [ cosd(obs_lat)*cosd(obs_lon);  cosd(obs_lat)*sind(obs_lon); sind(obs_lat)];

%% ── Pre-allocate outputs ─────────────────────────────────────────────────────
az       = zeros(1, N);
el       = zeros(1, N);
range_km = zeros(1, N);

%% ── Main loop over time steps ────────────────────────────────────────────────
for k = 1:N
    % Greenwich Mean Sidereal Time at this epoch (deg)
    GMST = mod(280.46061837 + 360.98564736629 * (jd(k) - 2451545.0), 360);
    cG = cosd(GMST);
    sG = sind(GMST);

    % Rotate satellite ECI -> ECEF: Rz(+GMST) * r_eci
    r_sat_k = r_sat_eci(:, k);
    r_sat_ecef = [ cG*r_sat_k(1) + sG*r_sat_k(2); ...
                  -sG*r_sat_k(1) + cG*r_sat_k(2); ...
                   r_sat_k(3)];

    % Range vector in ECEF
    rho_ecef = r_sat_ecef - r_obs_ecef;

    % Project onto ENU frame
    rho_E = dot(rho_ecef, e_E_ecef);
    rho_N = dot(rho_ecef, e_N_ecef);
    rho_U = dot(rho_ecef, e_U_ecef);

    % Range magnitude
    rng = norm(rho_ecef);

    % Elevation and azimuth
    el_k = asind(rho_U / rng);
    az_k = mod(atan2d(rho_E, rho_N), 360);

    az(k)       = az_k;
    el(k)       = el_k;
    range_km(k) = rng;
end
end
