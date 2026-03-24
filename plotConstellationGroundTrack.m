function fig = plotConstellationGroundTrack(sats, varargin)
%PLOTCONSTELLATIONGROUNDTRACK  Plot all constellation satellite ground tracks.
%
%   fig = plotConstellationGroundTrack(sats)
%   fig = plotConstellationGroundTrack(sats, Name, Value, ...)
%
%   Inputs:
%     sats - 1xT struct array of earthOrbit structs (from walkerConstellation)
%
%   Options (Name-Value pairs):
%     'NumOrbits'      - orbits to plot per satellite (default: 1)
%     'J2'             - apply J2 secular drift (default: true)
%     'ColorByPlane'   - color by plane_idx field if present (default: true)
%     'GroundStations' - struct array with .lat .lon .name (default: [])
%     'Title'          - override title string (default: auto-generated)

p = inputParser;
addParameter(p, 'NumOrbits',      1,    @isnumeric);
addParameter(p, 'J2',             true, @islogical);
addParameter(p, 'ColorByPlane',   true, @islogical);
addParameter(p, 'GroundStations', [],   @(x) isempty(x) || isstruct(x));
addParameter(p, 'Title',          '',   @ischar);
parse(p, varargin{:});
opts = p.Results;

mu_E = 398600.4418;
R_E  = 6378.1363;
J2c  = 1.08262668e-3;
om_E = 360.98564724 / 86400;   % deg/s (sidereal)

n_sats = numel(sats);

%% ── Determine plane assignments ─────────────────────────────────────────────
if opts.ColorByPlane && isfield(sats(1), 'plane_idx')
    plane_ids = [sats.plane_idx];
else
    % Group by RAAN value (round to nearest 0.01 deg for float safety)
    raan_vals    = round([sats.RAAN] * 100) / 100;
    unique_raans = unique(raan_vals);
    plane_ids    = zeros(1, n_sats);
    for k = 1:numel(unique_raans)
        plane_ids(raan_vals == unique_raans(k)) = k;
    end
end

P      = max(plane_ids);
cmap   = lines(P);

%% ── Figure setup (dark theme) ───────────────────────────────────────────────
bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];
axCol  = [0.14 0.18 0.26];

fig = figure('Color', bgCol, 'Position', [80 80 1300 620]);
ax  = axes('Parent', fig, 'Color', axCol, ...
    'XColor', txtCol * 0.7, 'YColor', txtCol * 0.7, ...
    'GridColor', [0.25 0.28 0.35], 'GridAlpha', 0.5);
hold(ax, 'on');  grid(ax, 'on');

%% ── Coastlines ───────────────────────────────────────────────────────────────
[coastlon, coastlat] = loadCoastlines();
if ~isempty(coastlon)
    plot(ax, coastlon, coastlat, '-', 'Color', [0.42 0.52 0.60], 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
end

%% ── Ground tracks for each satellite ────────────────────────────────────────
plane_plotted = false(1, P);   % track which planes have legend entries

for s = 1:n_sats
    orb_s  = sats(s);
    pid    = plane_ids(s);
    col    = cmap(pid, :);

    a_s     = orb_s.a;
    e_s     = orb_s.e;
    i_deg_s = orb_s.i;
    RAAN0_s = orb_s.RAAN;
    om0_s   = orb_s.omega;
    M0_s    = orb_s.M0;
    T_s     = orb_s.period;
    n_s     = 2*pi / T_s;

    % J2 secular drift rates (same as plotGroundTrack)
    if opts.J2
        p_orb    = a_s * (1 - e_s^2);
        factor   = 1.5 * n_s * J2c * (R_E / p_orb)^2;
        RAAN_dot = -factor * cosd(i_deg_s);
        om_dot   =  factor * (2.5 * cosd(i_deg_s)^2 - 1.0);
        n_corr   =  n_s + factor * sqrt(1-e_s^2) * (1.5 * cosd(i_deg_s)^2 - 0.5);
    else
        RAAN_dot = 0;  om_dot = 0;  n_corr = n_s;
    end

    % GMST at epoch
    GMST_epoch = mod(280.46061837 + 360.98564736629 * (orb_s.epoch_jd - 2451545.0), 360);

    % Time vector: opts.NumOrbits full periods, 360 steps per orbit
    N_steps = opts.NumOrbits * 360;
    t_vec   = linspace(0, opts.NumOrbits * T_s, N_steps + 1);

    lat_all = zeros(1, N_steps + 1);
    lon_all = zeros(1, N_steps + 1);

    for k = 1:N_steps+1
        t_k   = t_vec(k);
        RAAN  = RAAN0_s + rad2deg(RAAN_dot) * t_k;
        omega = om0_s   + rad2deg(om_dot)   * t_k;
        M_rad = mod(deg2rad(M0_s) + n_corr * t_k, 2*pi);
        E_k   = keplerSolve(M_rad, e_s);
        nu    = 2 * atan2(sqrt(1+e_s)*sin(E_k/2), sqrt(1-e_s)*cos(E_k/2));

        [r_eci, ~] = coe2eci(a_s, e_s, i_deg_s, RAAN, omega, rad2deg(nu));

        GMST_k = mod(GMST_epoch + om_E * t_k, 360);
        cG = cosd(-GMST_k);  sG = sind(-GMST_k);
        r_ecef = [cG*r_eci(1) - sG*r_eci(2);
                  sG*r_eci(1) + cG*r_eci(2);
                  r_eci(3)];

        lat_all(k) = asind(r_ecef(3) / norm(r_ecef));
        lon_all(k) = atan2d(r_ecef(2), r_ecef(1));
    end

    % Plot with antimeridian wrapping
    breaks = [0, find(abs(diff(lon_all)) > 180), numel(lon_all)];

    % Only show legend entry for first satellite in each plane
    show_legend = ~plane_plotted(pid);
    if show_legend
        plane_plotted(pid) = true;
        lbl = sprintf('Plane %d', pid);
        vis = 'on';
    else
        lbl = '';
        vis = 'off';
    end

    for seg_i = 1:numel(breaks)-1
        seg = (breaks(seg_i)+1) : breaks(seg_i+1);
        plot(ax, lon_all(seg), lat_all(seg), '-', ...
            'Color', [col, 0.7], ...
            'LineWidth', 0.8, ...
            'DisplayName', lbl, ...
            'HandleVisibility', iif(seg_i == 1 && show_legend, 'on', 'off'));
    end
end

%% ── Ground stations ──────────────────────────────────────────────────────────
if ~isempty(opts.GroundStations)
    for gs = opts.GroundStations(:)'
        scatter(ax, gs.lon, gs.lat, 100, 's', ...
            'MarkerFaceColor', [1.0 0.50 0.10], 'MarkerEdgeColor', 'none', ...
            'DisplayName', gs.name);
        text(ax, gs.lon + 2, gs.lat + 3, gs.name, ...
            'Color', [1.0 0.72 0.35], 'FontSize', 8);
    end
end

%% ── Axis formatting ──────────────────────────────────────────────────────────
xlim(ax, [-180 180]);  ylim(ax, [-90 90]);
xticks(ax, -180:30:180);  yticks(ax, -90:30:90);
xlabel(ax, 'Longitude (deg)', 'Color', txtCol, 'FontSize', 11);
ylabel(ax, 'Latitude (deg)',  'Color', txtCol, 'FontSize', 11);

% Title
if ~isempty(opts.Title)
    title_str = opts.Title;
else
    orb1 = sats(1);
    % Try to detect Walker info from constellation structure
    alt_km = round((orb1.a - R_E));

    % Check if we have plane_idx to detect T/P structure
    if isfield(sats(1), 'plane_idx')
        T_total = n_sats;
        P_total = max(plane_ids);
        % Try to infer F from first sat of second plane vs first plane
        if P_total > 1
            S = T_total / P_total;
            M_plane0_sat0 = sats(1).M0;
            % Find first sat in plane 2
            sat_p2 = find(plane_ids == 2, 1);
            if ~isempty(sat_p2)
                M_plane1_sat0 = sats(sat_p2).M0;
                % M_kj formula: M = (j*360/S + k*F*360/T) -> for j=0, k=1: M = F*360/T
                F_est = round(M_plane1_sat0 * T_total / 360);
                title_str = sprintf('Walker %d/%d/%d  |  %.0f km  |  i = %.1f°', ...
                    T_total, P_total, F_est, alt_km, orb1.i);
            else
                title_str = sprintf('Walker %d/%d  |  %.0f km  |  i = %.1f°', ...
                    T_total, P_total, alt_km, orb1.i);
            end
        else
            title_str = sprintf('Walker %d/1/0  |  %.0f km  |  i = %.1f°', ...
                n_sats, alt_km, orb1.i);
        end
    else
        title_str = sprintf('Constellation Ground Tracks (%d satellites)', n_sats);
    end
end

title(ax, title_str, 'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');

legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'TextColor', txtCol, 'Color', [0.08 0.10 0.14], 'EdgeColor', [0.30 0.32 0.38]);
end

%% ── Local: coastline loader (same pattern as plotGroundTrack) ────────────────
function [lon, lat] = loadCoastlines()
%LOADCOASTLINES  Load coastline data without requiring the Mapping Toolbox.
%   Tries (in order):
%     1. coastlines_cache.mat next to this file
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
    fprintf('plotConstellationGroundTrack: coastlines_cache.mat not found. Downloading now...\n');
    try
        downloadCoastlines();
        s = load(cachePath, 'coastlon', 'coastlat');
        lon = s.coastlon;  lat = s.coastlat;
    catch ME
        warning('plotConstellationGroundTrack:noCoastlines', ...
            'Could not load or download coastlines: %s\nPlotting without coastlines.', ME.message);
    end
end

%% ── Local helpers ────────────────────────────────────────────────────────────
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
