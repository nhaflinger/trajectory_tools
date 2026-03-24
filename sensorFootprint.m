function [fp, varargout] = sensorFootprint(lat_sat, lon_sat, alt_km, half_angle_deg, varargin)
%SENSORFOOTPRINT  Compute and optionally plot the ground footprint of a sensor.
%
%   fp = sensorFootprint(lat_sat, lon_sat, alt_km, half_angle_deg)
%   fp = sensorFootprint(lat_sat, lon_sat, alt_km, half_angle_deg, Name, Value, ...)
%   [fp, fig] = sensorFootprint(..., 'Plot', true)
%
%   Inputs:
%     lat_sat        - Sub-satellite latitude (deg, +North)
%     lon_sat        - Sub-satellite longitude (deg, +East)
%     alt_km         - Altitude above Earth surface (km)
%     half_angle_deg - Sensor half-cone angle from boresight direction (deg)
%
%   Name-Value options:
%     'PointAz_deg'  - Boresight azimuth offset (deg CW from North)  [default: 0]
%     'PointEl_deg'  - Off-nadir pointing angle (deg); 0 = nadir      [default: 0]
%     'N_pts'        - Number of boundary points                       [default: 360]
%     'Plot'         - Generate figure                                 [default: false]
%     'PlotAxes'     - Plot on existing axes handle                    [default: []]
%     'FillColor'    - RGB fill color for footprint polygon            [default: [0.30 0.75 0.93]]
%     'FillAlpha'    - Transparency of footprint fill                  [default: 0.3]
%     'ShowSubSat'   - Mark sub-satellite point                        [default: true]
%     'ShowCenter'   - Mark footprint center when off-nadir            [default: true]
%
%   Output struct fields:
%     lat_boundary          - Nx1 boundary latitudes (deg)
%     lon_boundary          - Nx1 boundary longitudes (deg)
%     lat_center            - Footprint center latitude (deg)
%     lon_center            - Footprint center longitude (deg)
%     lat_sat               - Sub-satellite latitude (deg)
%     lon_sat               - Sub-satellite longitude (deg)
%     alt_km                - Altitude (km)
%     half_angle_deg        - Sensor half-angle (deg)
%     earth_central_angle_deg - Earth central angle sub-sat to footprint edge (deg)
%     footprint_radius_km   - Approximate ground radius (km)
%     footprint_area_km2    - Approximate footprint area (km^2)
%     max_range_km          - Slant range to footprint edge (km)
%     min_elevation_deg     - Elevation angle at footprint edge (deg)

R_E = 6378.1363;   % km

%% ── Parse options ────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'PointAz_deg', 0,                    @isnumeric);
addParameter(p, 'PointEl_deg', 0,                    @isnumeric);
addParameter(p, 'N_pts',       360,                  @isnumeric);
addParameter(p, 'Plot',        false,                @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotAxes',    [],                   @(x) isempty(x) || ishandle(x));
addParameter(p, 'FillColor',   [0.30 0.75 0.93],     @isnumeric);
addParameter(p, 'FillAlpha',   0.3,                  @isnumeric);
addParameter(p, 'ShowSubSat',  true,                 @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ShowCenter',  true,                 @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
opts = p.Results;

doPlot = logical(opts.Plot) || ~isempty(opts.PlotAxes);

%% ── Satellite ECEF position ──────────────────────────────────────────────────
r_sat_mag = R_E + alt_km;
S_ecef    = r_sat_mag * [cosd(lat_sat)*cosd(lon_sat);
                          cosd(lat_sat)*sind(lon_sat);
                          sind(lat_sat)];   % 3x1

%% ── Nadir unit vector (pointing from satellite toward Earth center) ───────────
nadir_hat = -S_ecef / norm(S_ecef);   % 3x1

%% ── Local ENU basis at sub-satellite point ───────────────────────────────────
% East, North, Up unit vectors in ECEF
e_E = [-sind(lon_sat);             cosd(lon_sat);            0           ];
e_N = [-sind(lat_sat)*cosd(lon_sat); -sind(lat_sat)*sind(lon_sat); cosd(lat_sat)];
e_U = [ cosd(lat_sat)*cosd(lon_sat);  cosd(lat_sat)*sind(lon_sat); sind(lat_sat)];

%% ── Boresight direction ──────────────────────────────────────────────────────
pointEl = opts.PointEl_deg;
pointAz = opts.PointAz_deg;

if pointEl == 0
    % Pure nadir
    boresight_hat = nadir_hat;
else
    % Rotate nadir vector by pointEl toward the direction given by pointAz.
    % pointAz is clockwise from North in the local horizontal plane.
    % The horizontal direction in ECEF corresponding to pointAz:
    horiz_hat = cosd(pointAz) * e_N + sind(pointAz) * e_E;   % 3x1

    % Rotate nadir by pointEl around the axis perpendicular to (nadir, horiz)
    % Using Rodrigues' rotation: rotate nadir_hat by pointEl deg
    % toward horiz_hat (which is perpendicular to nadir_hat)
    boresight_hat = cosd(pointEl) * nadir_hat + sind(pointEl) * horiz_hat;
    boresight_hat = boresight_hat / norm(boresight_hat);
end

%% ── Footprint center (boresight intersection with Earth) ────────────────────
[hit_ok, P_center] = raySphereIntersect(S_ecef, boresight_hat, R_E);
if ~hit_ok
    error('sensorFootprint: boresight does not intersect Earth. Reduce PointEl_deg.');
end

lat_center = asind(P_center(3) / norm(P_center));
lon_center = atan2d(P_center(2), P_center(1));

%% ── Footprint boundary via ray-sphere intersection ───────────────────────────
N_pts     = opts.N_pts;
phi_vec   = linspace(0, 360*(1 - 1/N_pts), N_pts);   % azimuth around cone (deg)

lat_bnd = zeros(N_pts, 1);
lon_bnd = zeros(N_pts, 1);

% Build an orthonormal basis in the plane perpendicular to boresight_hat.
% Find two perpendicular vectors using Gram-Schmidt.
arb = [1; 0; 0];
if abs(dot(arb, boresight_hat)) > 0.9
    arb = [0; 1; 0];
end
perp1 = arb - dot(arb, boresight_hat) * boresight_hat;
perp1 = perp1 / norm(perp1);
perp2 = cross(boresight_hat, perp1);
perp2 = perp2 / norm(perp2);

for k = 1:N_pts
    phi = phi_vec(k);
    % Boundary ray direction: rotate boresight by half_angle_deg
    % around axis perp1*cos(phi) + perp2*sin(phi)
    % Using Rodrigues:  d = cos(ha)*bore + sin(ha)*(rot_axis x bore) ...
    % simpler: d = cos(ha)*bore + sin(ha)*(cos(phi)*perp1 + sin(phi)*perp2)
    d = cosd(half_angle_deg) * boresight_hat + ...
        sind(half_angle_deg) * (cosd(phi) * perp1 + sind(phi) * perp2);
    d = d / norm(d);

    [ok, P_hit] = raySphereIntersect(S_ecef, d, R_E);
    if ~ok
        % Ray missed Earth — clamp to Earth limb
        % Project d onto plane perpendicular to S_ecef
        d_tang = d - dot(d, e_U) * e_U;
        if norm(d_tang) < 1e-12
            d_tang = perp1;
        end
        d_tang = d_tang / norm(d_tang);
        % Horizon: rho = acos(R_E / r_sat_mag)
        rho_h = acosd(R_E / r_sat_mag);
        lat_bnd(k) = lat_sat + rho_h * (cosd(phi) * cosd(pointAz) - sind(phi) * sind(pointAz));
        lon_bnd(k) = lon_sat;
        warning('sensorFootprint:rayMiss', 'Boundary ray %d missed Earth. Clamping.', k);
        continue;
    end
    lat_bnd(k) = asind(P_hit(3) / norm(P_hit));
    lon_bnd(k) = atan2d(P_hit(2), P_hit(1));
end

% Close the polygon
lat_bnd(end+1) = lat_bnd(1);
lon_bnd(end+1) = lon_bnd(1);

%% ── Derived quantities (nadir-pointing reference geometry) ──────────────────
% Earth central angle for nadir pointing at half_angle_deg nadir angle
eta = half_angle_deg;  % nadir angle = sensor half-angle for nadir pointing
arg = (R_E + alt_km) / R_E * sind(eta);
if arg > 1
    rho_deg = 90 - eta;
else
    rho_deg = asind(arg) - eta;
end
earth_central_angle_deg = rho_deg;
footprint_radius_km     = R_E * deg2rad(rho_deg);

% Footprint area: spherical cap
footprint_area_km2 = 2 * pi * R_E^2 * (1 - cosd(rho_deg));

% Slant range and elevation at footprint edge
% Triangle: O, S, P  with angle at S = eta, OS = R_E+alt, OP = R_E
% rho = earth central angle; angle OPS = 90 + el (el = elevation at P)
% law of sines: sin(eta)/R_E = sin(rho)/range
range_edge_km   = R_E * sind(rho_deg) / sind(eta);
if eta < 1e-9
    range_edge_km = alt_km;
end
% Elevation at footprint edge
% el = pi - eta - rho - pi/2 = pi/2 - eta - rho... recompute properly
% angle at P = 180 - eta - rho (interior angles of triangle)
% elevation = angle_at_P - 90
el_edge_deg = 90 - eta - rho_deg;

%% ── Build output struct ──────────────────────────────────────────────────────
fp = struct( ...
    'lat_boundary',           lat_bnd,               ...
    'lon_boundary',           lon_bnd,               ...
    'lat_center',             lat_center,            ...
    'lon_center',             lon_center,            ...
    'lat_sat',                lat_sat,               ...
    'lon_sat',                lon_sat,               ...
    'alt_km',                 alt_km,                ...
    'half_angle_deg',         half_angle_deg,        ...
    'earth_central_angle_deg',earth_central_angle_deg, ...
    'footprint_radius_km',    footprint_radius_km,   ...
    'footprint_area_km2',     footprint_area_km2,    ...
    'max_range_km',           range_edge_km,         ...
    'min_elevation_deg',      el_edge_deg);

%% ── Plot ─────────────────────────────────────────────────────────────────────
if doPlot
    if ~isempty(opts.PlotAxes)
        ax  = opts.PlotAxes;
        fig = get(ax, 'Parent');
    else
        bgCol  = [0.10 0.12 0.16];
        txtCol = [0.88 0.88 0.88];
        axCol  = [0.14 0.18 0.26];

        fig = figure('Color', bgCol, 'Position', [100 80 1100 580]);
        ax  = axes('Parent', fig, 'Color', axCol, ...
            'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
            'GridColor', [0.25 0.28 0.35], 'GridAlpha', 0.5);
        hold(ax, 'on');  grid(ax, 'on');

        % Load and plot coastlines
        [clon, clat] = loadCoastlines();
        if ~isempty(clon)
            plot(ax, clon, clat, '-', 'Color', [0.42 0.52 0.60], ...
                'LineWidth', 0.7, 'HandleVisibility', 'off');
        end

        % Grid lines at 30 degrees
        for glon = -180:30:180
            xline(ax, glon, '-', 'Color', [0.20 0.24 0.32], 'LineWidth', 0.4, ...
                'HandleVisibility', 'off');
        end
        for glat = -90:30:90
            yline(ax, glat, '-', 'Color', [0.20 0.24 0.32], 'LineWidth', 0.4, ...
                'HandleVisibility', 'off');
        end

        xlim(ax, [-180 180]);  ylim(ax, [-90 90]);
        xticks(ax, -180:30:180);  yticks(ax, -90:30:90);
        xlabel(ax, 'Longitude (deg)', 'Color', txtCol, 'FontSize', 11);
        ylabel(ax, 'Latitude (deg)',  'Color', txtCol, 'FontSize', 11);

        ttl = sprintf('Sensor Footprint  |  Alt = %.0f km  |  Half-angle = %.1f°', ...
            alt_km, half_angle_deg);
        if pointEl > 0
            ttl = [ttl sprintf('  |  Off-nadir = %.1f° @ Az = %.0f°', pointEl, pointAz)];
        end
        title(ax, ttl, 'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');
    end

    hold(ax, 'on');

    % Handle antimeridian crossing: split polygon at ±180
    lon_b = lon_bnd;
    lat_b = lat_bnd;
    crosses = find(abs(diff(lon_b(1:end-1))) > 180);
    if isempty(crosses)
        fill(ax, lon_b, lat_b, opts.FillColor, ...
            'FaceAlpha', opts.FillAlpha, 'EdgeColor', opts.FillColor, ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(ax, lon_b, lat_b, '-', 'Color', opts.FillColor, ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    else
        % Split at each antimeridian crossing
        segs = [1; crosses+1; numel(lon_b)];
        for si = 1:numel(segs)-1
            idx = segs(si):segs(si+1);
            fill(ax, lon_b(idx), lat_b(idx), opts.FillColor, ...
                'FaceAlpha', opts.FillAlpha, 'EdgeColor', opts.FillColor, ...
                'LineWidth', 1.2, 'HandleVisibility', 'off');
            plot(ax, lon_b(idx), lat_b(idx), '-', 'Color', opts.FillColor, ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end

    % Sub-satellite point
    if logical(opts.ShowSubSat)
        scatter(ax, lon_sat, lat_sat, 80, 'x', ...
            'MarkerEdgeColor', [0.95 0.90 0.25], 'LineWidth', 2.0, ...
            'DisplayName', sprintf('Sub-sat (%.1f°N, %.1f°E)', lat_sat, lon_sat));
    end

    % Footprint center (only if different from sub-sat point)
    if logical(opts.ShowCenter) && pointEl > 0
        scatter(ax, lon_center, lat_center, 60, 'o', ...
            'MarkerFaceColor', opts.FillColor, 'MarkerEdgeColor', [0.95 0.95 0.95], ...
            'LineWidth', 1.0, ...
            'DisplayName', sprintf('FP center (%.1f°N, %.1f°E)', lat_center, lon_center));
    end

    % Annotation box — only when creating a new figure (not overlaying on existing axes)
    if isempty(opts.PlotAxes)
        bgCol_l  = [0.10 0.12 0.16];
        txtCol_l = [0.88 0.88 0.88];
        ann_str  = sprintf(['Alt: %.0f km\nHalf-angle: %.1f°\n' ...
            'FP radius: %.0f km\nMin el: %.1f°'], ...
            alt_km, half_angle_deg, footprint_radius_km, el_edge_deg);
        annotation(fig, 'textbox', [0.72 0.02 0.26 0.18], 'String', ann_str, ...
            'Color', txtCol_l, 'FontSize', 9, ...
            'BackgroundColor', bgCol_l, 'EdgeColor', [0.30 0.32 0.38]);
    end

    if nargout > 1
        varargout{1} = fig;
    end
else
    if nargout > 1
        varargout{1} = [];
    end
end

end  % main function

%% ── Local: ray–sphere intersection ─────────────────────────────────────────
function [ok, P] = raySphereIntersect(S, d, R)
%RAYSPHEREINTERSECT  Intersect ray S + t*d with sphere |P| = R.
%   Returns ok=true and the closer intersection point P if discriminant >= 0.
    Sd  = dot(S, d);
    disc = Sd^2 - (dot(S,S) - R^2);
    if disc < 0
        ok = false;  P = [nan; nan; nan];
        return;
    end
    t  = -Sd - sqrt(disc);   % smaller (closer) root
    if t < 0
        % Satellite is inside sphere (should not happen in normal usage)
        t = -Sd + sqrt(disc);
    end
    P  = S + t * d;
    ok = true;
end

%% ── Local: coastline loader ─────────────────────────────────────────────────
function [lon, lat] = loadCoastlines()
%LOADCOASTLINES  Load coastline data without requiring the Mapping Toolbox.
    lon = [];  lat = [];
    cachePath = fullfile(fileparts(mfilename('fullpath')), 'coastlines_cache.mat');
    if isfile(cachePath)
        s = load(cachePath, 'coastlon', 'coastlat');
        lon = s.coastlon;  lat = s.coastlat;
        return;
    end
    try
        s = load('coastlines');
        lon = s.coastlon;  lat = s.coastlat;
        return;
    catch
    end
    fprintf('sensorFootprint: coastlines_cache.mat not found. Downloading now...\n');
    try
        downloadCoastlines();
        s = load(cachePath, 'coastlon', 'coastlat');
        lon = s.coastlon;  lat = s.coastlat;
    catch ME
        warning('sensorFootprint:noCoastlines', ...
            'Could not load or download coastlines: %s\nPlotting without coastlines.', ME.message);
    end
end
