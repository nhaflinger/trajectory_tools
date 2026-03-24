function fig = plotGroundTrack(orb, varargin)
%PLOTGROUNDTRACK  Mercator ground track plot for an Earth orbit.
%
%   fig = plotGroundTrack(orb)
%   fig = plotGroundTrack(orb, Name, Value, ...)
%
%   orb must be a struct returned by earthOrbit().
%
%   Options (Name-Value pairs):
%     'NumOrbits'      - number of orbits to plot (default: 3)
%     'StepsPerOrbit'  - propagation steps per orbit (default: 360)
%     'J2'             - include J2 secular RAAN / argument-of-perigee drift (default: true)
%     'ColorByOrbit'   - color each orbit pass differently (default: true)
%     'GroundStations' - struct array with fields .lat, .lon, .name  (default: [])
%     'ShowNodes'      - mark ascending node crossings (default: true)
%
%   Returns: figure handle

p = inputParser;
addParameter(p, 'NumOrbits',     3,    @isnumeric);
addParameter(p, 'StepsPerOrbit', 360,  @isnumeric);
addParameter(p, 'J2',            true, @islogical);
addParameter(p, 'ColorByOrbit',  true, @islogical);
addParameter(p, 'GroundStations', [],  @(x) isempty(x) || isstruct(x));
addParameter(p, 'ShowNodes',     true, @islogical);
parse(p, varargin{:});
opts = p.Results;

mu_E = 398600.4418;
R_E  = 6378.1363;
J2   = 1.08262668e-3;
om_E = 360.98564724 / 86400;          % deg/s  (sidereal rotation rate)

a     = orb.a;
e     = orb.e;
i_deg = orb.i;
RAAN0 = orb.RAAN;
om0   = orb.omega;
M0    = orb.M0;
T     = orb.period;
n     = 2*pi / T;

% J2 secular drift rates (rad/s)
if opts.J2
    p_orb    = a * (1 - e^2);
    factor   = 1.5 * n * J2 * (R_E / p_orb)^2;
    RAAN_dot = -factor * cosd(i_deg);
    om_dot   =  factor * (2.5 * cosd(i_deg)^2 - 1.0);
    n_corr   =  n + factor * sqrt(1-e^2) * (1.5 * cosd(i_deg)^2 - 0.5);
else
    RAAN_dot = 0;  om_dot = 0;  n_corr = n;
end

% GMST at epoch (IAU approximation, deg)
GMST_epoch = mod(280.46061837 + 360.98564736629 * (orb.epoch_jd - 2451545.0), 360);

% Time vector
N_steps = opts.NumOrbits * opts.StepsPerOrbit;
t_vec   = linspace(0, opts.NumOrbits * T, N_steps + 1);

% Propagate: ECI -> ECEF -> geodetic (spherical)
lat_all   = zeros(1, N_steps + 1);
lon_all   = zeros(1, N_steps + 1);
orb_num   = zeros(1, N_steps + 1);

for k = 1:N_steps+1
    t     = t_vec(k);
    RAAN  = RAAN0 + rad2deg(RAAN_dot) * t;
    omega = om0   + rad2deg(om_dot)   * t;
    M     = mod(deg2rad(M0) + n_corr * t, 2*pi);
    E     = keplerSolve(M, e);
    nu    = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));

    [r_eci, ~] = coe2eci(a, e, i_deg, RAAN, omega, rad2deg(nu));

    % ECI -> ECEF: rotate around Z by -(GMST + om_E*t)
    GMST  = mod(GMST_epoch + om_E * t, 360);
    cG    = cosd(-GMST);  sG = sind(-GMST);
    r_ecef = [cG*r_eci(1) - sG*r_eci(2);
              sG*r_eci(1) + cG*r_eci(2);
              r_eci(3)];

    % ECEF -> geodetic (spherical; good to ~0.1 deg)
    lat_all(k) = asind(r_ecef(3) / norm(r_ecef));
    lon_all(k) = atan2d(r_ecef(2), r_ecef(1));
    orb_num(k) = floor(t / T) + 1;
end

%% ── Figure ─────────────────────────────────────────────────────────────────
bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];
axCol  = [0.14 0.18 0.26];

fig = figure('Color', bgCol, 'Position', [80 80 1300 620]);
ax  = axes('Parent', fig, 'Color', axCol, ...
    'XColor', txtCol * 0.7, 'YColor', txtCol * 0.7, ...
    'GridColor', [0.25 0.28 0.35], 'GridAlpha', 0.5);
hold(ax, 'on');  grid(ax, 'on');

%% ── Coastlines ──────────────────────────────────────────────────────────────
[coastlon, coastlat] = loadCoastlines();
if ~isempty(coastlon)
    plot(ax, coastlon, coastlat, '-', 'Color', [0.42 0.52 0.60], 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
end

%% ── Ground track (break line at antimeridian) ───────────────────────────────
N_orbs = opts.NumOrbits;
if opts.ColorByOrbit && N_orbs > 1
    cmap = parula(N_orbs);
else
    cmap = repmat([0.25 0.72 0.98], N_orbs, 1);
end

for on = 1:N_orbs
    idx = find(orb_num == on);
    if isempty(idx), continue; end
    lo = lon_all(idx);
    la = lat_all(idx);

    % Detect antimeridian crossings and split there
    breaks = [0, find(abs(diff(lo)) > 180), numel(lo)];
    col    = cmap(on, :);

    lbl = 'off';
    if on == 1, lbl = sprintf('Orbit 1-%d', N_orbs); end

    for s = 1:numel(breaks)-1
        seg = (breaks(s)+1) : breaks(s+1);
        plot(ax, lo(seg), la(seg), '-', 'Color', col, 'LineWidth', 1.3, ...
            'DisplayName', lbl, 'HandleVisibility', iif(s==1 && on==1, 'on', 'off'));
    end
end

%% ── Ascending node markers ──────────────────────────────────────────────────
if opts.ShowNodes
    first_node = true;
    for k = 2:N_steps+1
        if lat_all(k-1) < 0 && lat_all(k) >= 0
            dn = iif(first_node, 'Ascending node', 'off');
            scatter(ax, lon_all(k), lat_all(k), 30, '^', ...
                'MarkerEdgeColor', [0.95 0.90 0.25], ...
                'MarkerFaceColor', [0.95 0.90 0.25], ...
                'LineWidth', 0.5, 'DisplayName', dn, ...
                'HandleVisibility', iif(first_node,'on','off'));
            first_node = false;
        end
    end
end

%% ── Epoch position ──────────────────────────────────────────────────────────
scatter(ax, lon_all(1), lat_all(1), 90, 'o', ...
    'MarkerFaceColor', [0.20 1.00 0.45], 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Epoch position');

%% ── Ground stations ─────────────────────────────────────────────────────────
if ~isempty(opts.GroundStations)
    for gs = opts.GroundStations(:)'
        scatter(ax, gs.lon, gs.lat, 100, 's', ...
            'MarkerFaceColor', [1.0 0.50 0.10], 'MarkerEdgeColor', 'none', ...
            'DisplayName', gs.name);
        text(ax, gs.lon + 2, gs.lat + 3, gs.name, ...
            'Color', [1.0 0.72 0.35], 'FontSize', 8);
    end
end

%% ── Axis formatting ─────────────────────────────────────────────────────────
xlim(ax, [-180 180]);  ylim(ax, [-90 90]);
xticks(ax, -180:30:180);  yticks(ax, -90:30:90);
xlabel(ax, 'Longitude (deg)', 'Color', txtCol, 'FontSize', 11);
ylabel(ax, 'Latitude (deg)',  'Color', txtCol, 'FontSize', 11);

title_str = sprintf('%s | alt: %.0f – %.0f km | i = %.2f° | e = %.4f | T = %.1f min', ...
    upper(orb.type), orb.alt_peri, orb.alt_apo, orb.i, orb.e, orb.period/60);
title(ax, title_str, 'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');

legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'TextColor', txtCol, 'Color', [0.08 0.10 0.14], 'EdgeColor', [0.30 0.32 0.38]);
end

%% ── Local helper: coastline loader ─────────────────────────────────────────
function [lon, lat] = loadCoastlines()
%LOADCOASTLINES  Load coastline data without requiring the Mapping Toolbox.
%   Tries (in order):
%     1. coastlines_cache.mat next to plotGroundTrack.m
%     2. Mapping Toolbox built-in 'coastlines'
%     3. Auto-download from Natural Earth via downloadCoastlines()
    lon = [];  lat = [];
    cachePath = fullfile(fileparts(mfilename('fullpath')), 'coastlines_cache.mat');
    if isfile(cachePath)
        s = load(cachePath, 'coastlon', 'coastlat');
        lon = s.coastlon;  lat = s.coastlat;
        return;
    end
    try
        s = load('coastlines');   % Mapping Toolbox
        lon = s.coastlon;  lat = s.coastlat;
        return;
    catch
    end
    % Neither source found — try to download automatically
    fprintf('plotGroundTrack: coastlines_cache.mat not found. Downloading now...\n');
    try
        downloadCoastlines();
        s = load(cachePath, 'coastlon', 'coastlat');
        lon = s.coastlon;  lat = s.coastlat;
    catch ME
        warning('plotGroundTrack:noCoastlines', ...
            'Could not load or download coastlines: %s\nPlotting without coastlines.', ME.message);
    end
end

%% ── Local helper: ternary ───────────────────────────────────────────────────
function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end

function E = keplerSolve(M, e)
    E = M;
    for k = 1:50
        dE = (M - E + e*sin(E)) / (1 - e*cos(E));
        E  = E + dE;
        if abs(dE) < 1e-13, break; end
    end
end
