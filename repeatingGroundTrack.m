function result = repeatingGroundTrack(N_rev, N_day, inc_deg, varargin)
%REPEATINGGROUNDTRACK  Find semi-major axis for a repeating ground track orbit.
%
%   result = repeatingGroundTrack(N_rev, N_day, inc_deg)
%   result = repeatingGroundTrack(N_rev, N_day, inc_deg, Name, Value, ...)
%
%   The satellite completes exactly N_rev revolutions in N_day sidereal days,
%   then repeats its ground track.
%
%   Inputs:
%     N_rev   - number of satellite revolutions in repeat cycle
%     N_day   - number of sidereal days in repeat cycle
%     inc_deg - inclination (deg)
%
%   Options (Name-Value pairs):
%     'e' - eccentricity (default: 0)
%
%   Returns struct with fields:
%     a_km, alt_km, e, i_deg, N_rev, N_day, period_s,
%     RAAN_dot_deg_day, ground_track_spacing_deg, orb

p = inputParser;
addParameter(p, 'e', 0, @isnumeric);
parse(p, varargin{:});
ecc = p.Results.e;

%% ── Constants ───────────────────────────────────────────────────────────────
mu_E       = 398600.4418;       % km^3/s^2
R_E        = 6378.1363;         % km
J2         = 1.08262668e-3;
omega_E    = 7.2921150e-5;      % rad/s
T_sidereal = 86164.1;           % s

%% ── Iterative J2-corrected solution ────────────────────────────────────────
% Initial Keplerian guess (no J2):
%   N_rev * n = N_day * omega_E
%   n = N_day * omega_E / N_rev
%   a = (mu / n^2)^(1/3)
n_init = N_day * omega_E / N_rev;   % but use sidereal day properly:
% Repeat condition: N_rev * n_M = N_day * (omega_E - RAAN_dot)
% For Keplerian start: assume RAAN_dot = 0
a = (mu_E * (N_day * T_sidereal / (N_rev * 2*pi))^2)^(1/3);

for iter = 1:20
    a_prev = a;
    n_M    = sqrt(mu_E / a^3);   % mean motion (rad/s)

    % J2 RAAN drift
    RAAN_dot = -1.5 * n_M * J2 * (R_E / a)^2 * cosd(inc_deg) / (1 - ecc^2)^2;

    % Required mean motion to satisfy repeat condition
    n_target = N_day * (omega_E - RAAN_dot) / N_rev;

    % New semi-major axis
    a = (mu_E / n_target^2)^(1/3);

    if abs(a - a_prev) < 1e-6
        break;
    end
end

%% ── Derived quantities ──────────────────────────────────────────────────────
R_E_val = 6378.1363;
alt_km  = a - R_E_val;
period_s = 2*pi * sqrt(a^3 / mu_E);

n_M_final   = sqrt(mu_E / a^3);
RAAN_dot_final = -1.5 * n_M_final * J2 * (R_E / a)^2 * cosd(inc_deg) / (1 - ecc^2)^2;
RAAN_dot_deg_day = rad2deg(RAAN_dot_final) * 86400;

ground_track_spacing_deg = 360 * N_day / N_rev;

% Create earthOrbit struct
orb = earthOrbit('coe', a, ecc, inc_deg, 0, 0, 0);

%% ── Output struct ────────────────────────────────────────────────────────────
result = struct( ...
    'a_km',                    a,                         ...
    'alt_km',                  alt_km,                    ...
    'e',                       ecc,                       ...
    'i_deg',                   inc_deg,                   ...
    'N_rev',                   N_rev,                     ...
    'N_day',                   N_day,                     ...
    'period_s',                period_s,                  ...
    'RAAN_dot_deg_day',        RAAN_dot_deg_day,          ...
    'ground_track_spacing_deg',ground_track_spacing_deg,  ...
    'orb',                     orb);

%% ── Print summary ───────────────────────────────────────────────────────────
fprintf('\n=== Repeating Ground Track: %d rev / %d day ===\n', N_rev, N_day);
fprintf('  Repeat ratio     : %d / %d (%.4f rev/day)\n', N_rev, N_day, N_rev/N_day);
fprintf('  Inclination      : %.2f deg\n',   inc_deg);
fprintf('  Semi-major axis  : %.4f km\n',    a);
fprintf('  Altitude         : %.4f km\n',    alt_km);
fprintf('  Period           : %.4f min\n',   period_s/60);
fprintf('  RAAN drift       : %.6f deg/day\n', RAAN_dot_deg_day);
fprintf('  Track spacing    : %.4f deg\n',   ground_track_spacing_deg);
fprintf('\n');
end
