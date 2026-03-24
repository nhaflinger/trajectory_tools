function res = hohmannTransfer(arg1, arg2, varargin)
%HOHMANNTRANSFER  Hohmann (or bi-elliptic) transfer between two circular Earth orbits.
%
%   res = hohmannTransfer(alt1_km, alt2_km)
%   res = hohmannTransfer(orb1, alt2_km)
%   res = hohmannTransfer(orb1, orb2)
%   res = hohmannTransfer(..., 'Bielliptic', r_int_km)
%
%   Each positional argument can be either:
%     - A scalar altitude in km above Earth surface, OR
%     - An earthOrbit struct (uses orb.a to extract semi-major axis)
%
%   Both orbits are treated as circular.
%
%   Name-Value options:
%     'Bielliptic'  - intermediate apoapsis radius in km (measured from Earth center).
%                     When provided, a bi-elliptic transfer is computed instead.
%
%   Output struct fields:
%     dv1_km_s, dv2_km_s, dv_total_km_s
%     tof_s, tof_min
%     alt1_km, alt2_km, r1_km, r2_km
%     a_transfer_km   (Hohmann) or a1_km, a2_km (bi-elliptic)
%     type            'hohmann' or 'bielliptic'
%     (bi-elliptic only) dv3_km_s

mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km

%% ── Parse positional arguments ─────────────────────────────────────────────
r1 = extractRadius(arg1, R_E);
r2 = extractRadius(arg2, R_E);

alt1 = r1 - R_E;
alt2 = r2 - R_E;

%% ── Parse name-value options ────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'Bielliptic', [], @(x) isempty(x) || isnumeric(x));
parse(p, varargin{:});
r_int = p.Results.Bielliptic;

%% ── Compute transfer ────────────────────────────────────────────────────────
if isempty(r_int)
    %% Hohmann transfer
    a_t = (r1 + r2) / 2;

    v_circ1   = sqrt(mu_E / r1);
    v_circ2   = sqrt(mu_E / r2);
    v_t_at_r1 = sqrt(mu_E * (2/r1 - 1/a_t));
    v_t_at_r2 = sqrt(mu_E * (2/r2 - 1/a_t));

    dv1 = v_t_at_r1 - v_circ1;   % negative means retrograde
    dv2 = v_circ2   - v_t_at_r2;  % negative means retrograde

    dv_total = abs(dv1) + abs(dv2);
    tof      = pi * sqrt(a_t^3 / mu_E);

    res = struct( ...
        'type',          'hohmann',   ...
        'alt1_km',       alt1,        ...
        'alt2_km',       alt2,        ...
        'r1_km',         r1,          ...
        'r2_km',         r2,          ...
        'a_transfer_km', a_t,         ...
        'dv1_km_s',      dv1,         ...
        'dv2_km_s',      dv2,         ...
        'dv_total_km_s', dv_total,    ...
        'tof_s',         tof,         ...
        'tof_min',       tof/60);

    fprintf('\n=== Hohmann Transfer ===\n');
    fprintf('  Orbit 1 altitude : %10.2f km  (r = %.2f km)\n', alt1, r1);
    fprintf('  Orbit 2 altitude : %10.2f km  (r = %.2f km)\n', alt2, r2);
    fprintf('  Transfer SMA     : %10.2f km\n', a_t);
    fprintf('  DV1              : %+10.4f km/s  (%+.2f m/s)\n', dv1, dv1*1000);
    fprintf('  DV2              : %+10.4f km/s  (%+.2f m/s)\n', dv2, dv2*1000);
    fprintf('  Total |DV|       : %10.4f km/s  (%.2f m/s)\n', dv_total, dv_total*1000);
    fprintf('  TOF              : %10.2f s  (%.2f min)\n', tof, tof/60);

else
    %% Bi-elliptic transfer
    if r_int < max(r1, r2)
        error('hohmannTransfer: bi-elliptic intermediate radius (%.1f km) must be >= max(r1,r2) = %.1f km', ...
            r_int, max(r1,r2));
    end

    a1 = (r1    + r_int) / 2;
    a2 = (r2    + r_int) / 2;

    dv1 = sqrt(mu_E*(2/r1    - 1/a1)) - sqrt(mu_E/r1);
    dv2 = sqrt(mu_E*(2/r_int - 1/a2)) - sqrt(mu_E*(2/r_int - 1/a1));
    dv3 = sqrt(mu_E/r2)                - sqrt(mu_E*(2/r2    - 1/a2));

    dv_total = abs(dv1) + abs(dv2) + abs(dv3);
    tof      = pi*(sqrt(a1^3/mu_E) + sqrt(a2^3/mu_E));

    res = struct( ...
        'type',          'bielliptic', ...
        'alt1_km',       alt1,         ...
        'alt2_km',       alt2,         ...
        'r1_km',         r1,           ...
        'r2_km',         r2,           ...
        'a1_km',         a1,           ...
        'a2_km',         a2,           ...
        'dv1_km_s',      dv1,          ...
        'dv2_km_s',      dv2,          ...
        'dv3_km_s',      dv3,          ...
        'dv_total_km_s', dv_total,     ...
        'tof_s',         tof,          ...
        'tof_min',       tof/60);

    fprintf('\n=== Bi-elliptic Transfer ===\n');
    fprintf('  Orbit 1 altitude   : %10.2f km  (r = %.2f km)\n', alt1, r1);
    fprintf('  Orbit 2 altitude   : %10.2f km  (r = %.2f km)\n', alt2, r2);
    fprintf('  Intermediate r_int : %10.2f km\n', r_int);
    fprintf('  Transfer ellipse 1 SMA : %10.2f km\n', a1);
    fprintf('  Transfer ellipse 2 SMA : %10.2f km\n', a2);
    fprintf('  DV1                : %+10.4f km/s  (%+.2f m/s)\n', dv1, dv1*1000);
    fprintf('  DV2                : %+10.4f km/s  (%+.2f m/s)\n', dv2, dv2*1000);
    fprintf('  DV3                : %+10.4f km/s  (%+.2f m/s)\n', dv3, dv3*1000);
    fprintf('  Total |DV|         : %10.4f km/s  (%.2f m/s)\n', dv_total, dv_total*1000);
    fprintf('  TOF                : %10.2f s  (%.2f h)\n', tof, tof/3600);
end

end

%% ── Local helper ─────────────────────────────────────────────────────────────
function r = extractRadius(arg, R_E)
% Return geocentric radius (km) from either a scalar altitude or an earthOrbit struct.
    if isstruct(arg)
        % earthOrbit struct: use semi-major axis as radius (circular approximation)
        r = arg.a;
    else
        r = R_E + arg;
    end
end
