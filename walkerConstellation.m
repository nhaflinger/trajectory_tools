function [sats, info] = walkerConstellation(inc_deg, alt_km, T, P, F, varargin)
%WALKERCONSTELLATION  Generate a Walker delta constellation (T/P/F notation).
%
%   [sats, info] = walkerConstellation(inc_deg, alt_km, T, P, F)
%   [sats, info] = walkerConstellation(inc_deg, alt_km, T, P, F, Name, Value, ...)
%
%   Inputs:
%     inc_deg - inclination (deg)
%     alt_km  - circular orbit altitude (km)
%     T       - total number of satellites
%     P       - number of orbital planes (must divide T evenly)
%     F       - phasing parameter (integer, 0 <= F <= P-1)
%
%   Options (Name-Value pairs):
%     'epoch_jd' - epoch for all orbits (default: 2451545.0, J2000)
%
%   Returns:
%     sats - 1xT struct array of earthOrbit structs with extra fields:
%              plane_idx, sat_idx_in_plane, sat_idx_total
%     info - struct with constellation parameters

p = inputParser;
addParameter(p, 'epoch_jd', 2451545.0, @isnumeric);
parse(p, varargin{:});
opts = p.Results;

%% ── Input validation ────────────────────────────────────────────────────────
if mod(T, P) ~= 0
    error('walkerConstellation: T (%d) must be divisible by P (%d)', T, P);
end
if F < 0 || F > P-1
    error('walkerConstellation: F (%d) must satisfy 0 <= F <= P-1 (%d)', F, P-1);
end

mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km

S    = T / P;         % satellites per plane
a_km = R_E + alt_km;  % semi-major axis

period_s = 2*pi * sqrt(a_km^3 / mu_E);

%% ── Build satellite array ───────────────────────────────────────────────────
% Pre-allocate using a template that already has the Walker-specific fields.
% repmat requires all elements to have identical fields, so the extra fields
% must be present in the template before the array is created.
orb_template = earthOrbit('coe', a_km, 0, inc_deg, 0, 0, 0, opts.epoch_jd);
orb_template.plane_idx        = 0;
orb_template.sat_idx_in_plane = 0;
orb_template.sat_idx_total    = 0;
sats = repmat(orb_template, 1, T);

sat_total = 0;
for k = 0:P-1          % plane index (0-based)
    RAAN_k = k * 360 / P;

    for j = 0:S-1      % satellite index within plane (0-based)
        sat_total = sat_total + 1;
        M_kj = mod(j * 360/S + k * F * 360/T, 360);

        orb_kj = earthOrbit('coe', a_km, 0, inc_deg, RAAN_k, 0, M_kj, opts.epoch_jd);

        % Append Walker-specific fields
        orb_kj.plane_idx         = k + 1;          % 1-based
        orb_kj.sat_idx_in_plane  = j + 1;          % 1-based
        orb_kj.sat_idx_total     = sat_total;       % 1-based

        sats(sat_total) = orb_kj;
    end
end

%% ── Info struct ─────────────────────────────────────────────────────────────
notation = sprintf('%d/%d/%d', T, P, F);

info = struct( ...
    'T',        T,              ...
    'P',        P,              ...
    'F',        F,              ...
    'S',        S,              ...
    'inc_deg',  inc_deg,        ...
    'alt_km',   alt_km,         ...
    'a_km',     a_km,           ...
    'period_s', period_s,       ...
    'notation', notation);

%% ── Print summary ───────────────────────────────────────────────────────────
fprintf('\n=== Walker Constellation %s ===\n', notation);
fprintf('  Altitude         : %.0f km\n',    alt_km);
fprintf('  Inclination      : %.2f deg\n',   inc_deg);
fprintf('  Total satellites : %d\n',          T);
fprintf('  Orbital planes   : %d\n',          P);
fprintf('  Sats per plane   : %d\n',          S);
fprintf('  Phasing param F  : %d\n',          F);
fprintf('  RAAN spacing     : %.2f deg\n',    360/P);
fprintf('  Period           : %.2f min\n',    period_s/60);
fprintf('  Semi-major axis  : %.2f km\n',     a_km);
fprintf('\n');
end
