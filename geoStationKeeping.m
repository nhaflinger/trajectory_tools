function budget = geoStationKeeping(varargin)
%GEOSTATIONKEEPING  Compute annual station-keeping delta-V budget for a GEO satellite.
%
%   budget = geoStationKeeping()
%   budget = geoStationKeeping(lon_deg)
%   budget = geoStationKeeping(lon_deg, year)
%
%   Inputs:
%     lon_deg  - geostationary longitude slot (deg East, default: 0)
%     year     - year for reference (default: current year)
%
%   North-South station keeping:
%     Lunar and solar gravity cause the orbital inclination to grow at
%     ~0.75 deg/year. The N-S DV budget dominates for most GEO satellites.
%
%   East-West station keeping:
%     Earth's triaxiality causes longitude drift toward stable equilibria
%     near 75 deg E and 255 deg E. The drift acceleration depends on
%     the satellite's longitude slot.
%
%   Output struct fields:
%     dv_ns_km_s           - N-S delta-V per year (km/s)
%     dv_ns_m_s            - N-S delta-V per year (m/s)
%     dv_ew_km_s           - E-W delta-V per year (km/s)
%     dv_ew_m_s            - E-W delta-V per year (m/s)
%     dv_total_km_s        - total delta-V per year (km/s)
%     dv_total_m_s         - total delta-V per year (m/s)
%     lon_deg              - longitude slot (deg E)
%     v_geo_km_s           - GEO circular velocity (km/s)
%     a_geo_km             - GEO semi-major axis (km)
%     lifetime_yr_per_10ms - years of station keeping per 10 m/s of propellant
%     deadband_ew_deg      - assumed E-W deadband (deg)

%% ── Constants ────────────────────────────────────────────────────────────────
mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km
om_E = 7.2921150e-5;  % rad/s

%% ── Parse inputs ─────────────────────────────────────────────────────────────
lon_deg = 0;
year    = clock;
year    = year(1);

if nargin >= 1, lon_deg = varargin{1}; end
if nargin >= 2, year    = varargin{2}; end

%% ── GEO orbit parameters ─────────────────────────────────────────────────────
a_geo     = (mu_E / om_E^2)^(1/3);    % ~42164.17 km
v_geo     = sqrt(mu_E / a_geo);        % ~3.0747 km/s

%% ── North-South station keeping (inclination control) ────────────────────────
% Lunar-solar perturbations cause inclination to grow at ~0.75 deg/year (mean).
% This is the dominant SK cost for most GEO satellites.
di_per_year_deg = 0.75;   % deg/year (representative mean over 18.6-yr cycle)

% DV_NS = v_GEO * delta_i_rad  (for small angles, 2*v*sin(di/2) ≈ v*di)
dv_ns_km_s = v_geo * deg2rad(di_per_year_deg);   % km/s/year

%% ── East-West station keeping (longitude control) ────────────────────────────
% Earth's tesseral harmonics (triaxiality) produce a longitude-dependent
% drift acceleration. Stable equilibria near 75 deg E and 255 deg E.
% Drift acceleration (deg/day^2):
%   a_lambda = -0.00168 * sin(2*(lon - 75))
a_lambda_deg_day2 = -0.00168 * sind(2 * (lon_deg - 75));   % deg/day^2

% Convert to km/day^2 at GEO radius: delta_v = a * R_GEO [km] in arc units
% More precisely: tangential acceleration = a_lambda [deg/day^2] * (pi/180) * a_geo [km/day^2 in arc length]
% But we need km/s^2. Convert:
%   a_km_day2 = a_lambda_deg_day2 * (pi/180) * a_geo   [km/day^2]
a_km_day2 = abs(a_lambda_deg_day2) * (pi/180) * a_geo;   % km/day^2

% Typical E-W deadband
deadband_ew_deg = 0.05;   % deg, typical commercial GEO

% E-W DV calculation using maneuver cycle approach:
% If drift acceleration is a [km/day^2], deadband is d [km] at GEO radius:
%   d_km = deadband_ew_deg * (pi/180) * a_geo
%   time to drift across deadband: t = sqrt(2*d_km / a_km_day2)  [days]
%   DV per cycle: dv_cycle = a_km_day2 * t  [km/day] -> convert to km/s
%   Cycles per year: N = T_year / t
%   DV_EW per year = N * dv_cycle
% This simplifies to: DV_EW = 2 * sqrt(2 * a_km_day2 * d_km) / 86400  [km/s/year]
d_km = deadband_ew_deg * (pi/180) * a_geo;   % km of arc at GEO radius

if a_km_day2 < 1e-12
    % At or very near a stable equilibrium point: negligible E-W SK cost
    dv_ew_km_s = 0;
else
    T_year_days = 365.25;
    % t_cycle = sqrt(2 * d_km / a_km_day2)   [days]
    t_cycle_days = sqrt(2 * d_km / a_km_day2);   % days
    dv_cycle_km_day = a_km_day2 * t_cycle_days;  % km/day
    dv_cycle_km_s   = dv_cycle_km_day / 86400;   % km/s (1 day = 86400 s)
    n_cycles = T_year_days / t_cycle_days;
    dv_ew_km_s = 2 * n_cycles * dv_cycle_km_s;   % factor 2: correct each cycle
end

%% ── Total budget ─────────────────────────────────────────────────────────────
dv_total_km_s = dv_ns_km_s + dv_ew_km_s;

% Lifetime sizing metric: years of SK life per 10 m/s of propellant budget
lifetime_yr_per_10ms = (10e-3) / dv_total_km_s;   % 10 m/s = 0.010 km/s

%% ── Print summary ────────────────────────────────────────────────────────────
fprintf('\n=== GEO Station-Keeping Budget ===\n');
fprintf('  Longitude slot    : %.1f deg E\n', lon_deg);
fprintf('  GEO radius        : %.2f km\n', a_geo);
fprintf('  GEO velocity      : %.4f km/s\n', v_geo);
fprintf('  Year              : %d\n', round(year));
fprintf('\n');
fprintf('  %-35s  %8s  %8s\n', 'Component', 'km/s/yr', 'm/s/yr');
fprintf('  %s\n', repmat('-', 1, 57));
fprintf('  %-35s  %8.4f  %8.2f\n', 'N-S (inclination, lunar-solar)', dv_ns_km_s, dv_ns_km_s*1e3);
fprintf('  %-35s  %8.4f  %8.2f\n', 'E-W (longitude, triaxiality)', dv_ew_km_s, dv_ew_km_s*1e3);
fprintf('  %s\n', repmat('-', 1, 57));
fprintf('  %-35s  %8.4f  %8.2f\n', 'TOTAL', dv_total_km_s, dv_total_km_s*1e3);
fprintf('\n');
fprintf('  E-W deadband assumed  : +/-%.3f deg\n', deadband_ew_deg);
fprintf('  E-W drift accel       : %.6f deg/day^2\n', a_lambda_deg_day2);
fprintf('  Lifetime per 10 m/s   : %.2f years\n', lifetime_yr_per_10ms);
fprintf('\n');

%% ── Build output struct ──────────────────────────────────────────────────────
budget = struct( ...
    'dv_ns_km_s',           dv_ns_km_s,           ...
    'dv_ns_m_s',            dv_ns_km_s * 1e3,     ...
    'dv_ew_km_s',           dv_ew_km_s,           ...
    'dv_ew_m_s',            dv_ew_km_s * 1e3,     ...
    'dv_total_km_s',        dv_total_km_s,        ...
    'dv_total_m_s',         dv_total_km_s * 1e3,  ...
    'lon_deg',              lon_deg,              ...
    'v_geo_km_s',           v_geo,                ...
    'a_geo_km',             a_geo,                ...
    'lifetime_yr_per_10ms', lifetime_yr_per_10ms, ...
    'deadband_ew_deg',      deadband_ew_deg,      ...
    'a_lambda_deg_day2',    a_lambda_deg_day2,    ...
    'di_per_year_deg',      di_per_year_deg);
end
