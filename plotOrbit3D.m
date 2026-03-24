function fig = plotOrbit3D(orb, varargin)
%PLOTORBIT3D  Three-dimensional ECI visualization of an Earth orbit.
%
%   fig = plotOrbit3D(orb)
%   fig = plotOrbit3D(orb, Name, Value, ...)
%
%   orb must be a struct returned by earthOrbit().
%
%   Options (Name-Value pairs):
%     'NumOrbits'     - number of orbit revolutions to trace (default: 1)
%     'StepsPerOrbit' - propagation resolution (default: 720)
%     'J2'            - include J2 secular drift (default: false)
%     'ShowEquator'   - translucent equatorial plane (default: true)
%     'ShowNodeLine'  - dashed line through ascending / descending nodes (default: true)
%     'ShowPerigee'   - mark perigee (for e > 0.01) (default: true)
%
%   Returns: figure handle

p = inputParser;
addParameter(p, 'NumOrbits',     1,     @isnumeric);
addParameter(p, 'StepsPerOrbit', 720,   @isnumeric);
addParameter(p, 'J2',            false, @islogical);
addParameter(p, 'ShowEquator',   true,  @islogical);
addParameter(p, 'ShowNodeLine',  true,  @islogical);
addParameter(p, 'ShowPerigee',   true,  @islogical);
parse(p, varargin{:});
opts = p.Results;

mu_E = 398600.4418;
R_E  = 6378.1363;
J2   = 1.08262668e-3;

a     = orb.a;
e     = orb.e;
i_deg = orb.i;
RAAN0 = orb.RAAN;
om0   = orb.omega;
M0    = orb.M0;
T     = orb.period;
n     = 2*pi / T;

% J2 secular drift rates
if opts.J2
    p_orb    = a * (1 - e^2);
    factor   = 1.5 * n * J2 * (R_E / p_orb)^2;
    RAAN_dot = -factor * cosd(i_deg);
    om_dot   =  factor * (2.5 * cosd(i_deg)^2 - 1.0);
    n_corr   =  n + factor * sqrt(1-e^2) * (1.5 * cosd(i_deg)^2 - 0.5);
else
    RAAN_dot = 0;  om_dot = 0;  n_corr = n;
end

% Propagate
N_steps = opts.NumOrbits * opts.StepsPerOrbit;
t_vec   = linspace(0, opts.NumOrbits * T, N_steps + 1);
xyz     = zeros(3, N_steps + 1);

for k = 1:N_steps+1
    t     = t_vec(k);
    RAAN  = RAAN0 + rad2deg(RAAN_dot) * t;
    omega = om0   + rad2deg(om_dot)   * t;
    M     = mod(deg2rad(M0) + n_corr * t, 2*pi);
    E     = keplerSolve(M, e);
    nu    = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    [r_eci, ~] = coe2eci(a, e, i_deg, RAAN, omega, rad2deg(nu));
    xyz(:, k) = r_eci;
end

r_max = max(vecnorm(xyz));

%% ── Figure ──────────────────────────────────────────────────────────────────
bgCol = [0.07 0.09 0.13];
fig = figure('Color', bgCol, 'Position', [120 80 920 820]);
ax  = axes('Parent', fig, 'Color', bgCol, ...
    'XColor', [0.55 0.55 0.60], 'YColor', [0.55 0.55 0.60], 'ZColor', [0.55 0.55 0.60]);
hold(ax, 'on');  grid(ax, 'on');  axis(ax, 'equal');
view(ax, 35, 22);
ax.GridColor = [0.20 0.22 0.28];  ax.GridAlpha = 0.5;

%% ── Earth sphere ─────────────────────────────────────────────────────────────
[xs, ys, zs] = sphere(72);
surf(ax, R_E*xs, R_E*ys, R_E*zs, ...
    'FaceColor', [0.12 0.28 0.50], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.90, 'DisplayName', 'Earth');

% Equatorial ring on the surface
th = linspace(0, 2*pi, 361);
plot3(ax, R_E*cos(th), R_E*sin(th), zeros(size(th)), ...
    '-', 'Color', [0.40 0.60 0.90 0.55], 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');

% Prime meridian line (0° longitude reference)
phi = linspace(-pi/2, pi/2, 181);
plot3(ax, R_E*cos(phi), zeros(size(phi)), R_E*sin(phi), ...
    ':', 'Color', [0.40 0.60 0.90 0.40], 'LineWidth', 0.6, ...
    'HandleVisibility', 'off');

%% ── Equatorial plane ─────────────────────────────────────────────────────────
if opts.ShowEquator
    r_pl = r_max * 1.08;
    [xe, ye] = meshgrid([-r_pl r_pl]);
    surf(ax, xe, ye, zeros(2,2), ...
        'FaceColor', [0.18 0.38 0.70], 'FaceAlpha', 0.07, 'EdgeColor', 'none', ...
        'DisplayName', 'Equatorial plane');
    % Ring at orbit apogee distance for scale reference
    r_ref = max(a*(1+e), R_E) * 1.0;
    plot3(ax, r_ref*cos(th), r_ref*sin(th), zeros(size(th)), ...
        '-', 'Color', [0.25 0.40 0.65 0.30], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
end

%% ── Node line ────────────────────────────────────────────────────────────────
if opts.ShowNodeLine
    r_nd  = r_max * 1.12;
    RAAN_r = deg2rad(RAAN0);
    plot3(ax, r_nd*[-cos(RAAN_r), cos(RAAN_r)], ...
              r_nd*[-sin(RAAN_r), sin(RAAN_r)], [0 0], ...
        '--', 'Color', [0.80 0.80 0.35 0.65], 'LineWidth', 1.0, ...
        'DisplayName', 'Line of nodes');
end

%% ── Orbit trace (colored by time: blue -> yellow) ────────────────────────────
N    = N_steps + 1;
cmap = parula(N);
for k = 1:N-1
    plot3(ax, xyz(1,k:k+1), xyz(2,k:k+1), xyz(3,k:k+1), '-', ...
        'Color', cmap(k,:), 'LineWidth', 1.6, 'HandleVisibility', 'off');
end
% Legend proxy
plot3(ax, NaN, NaN, NaN, '-', 'Color', cmap(round(N/2),:), 'LineWidth', 1.6, ...
    'DisplayName', 'Orbit trace');

%% ── Epoch position ───────────────────────────────────────────────────────────
scatter3(ax, orb.r_vec(1), orb.r_vec(2), orb.r_vec(3), 90, 'o', ...
    'MarkerFaceColor', [0.20 1.00 0.45], 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Epoch position');

%% ── Perigee marker ───────────────────────────────────────────────────────────
if opts.ShowPerigee && e > 0.01
    [r_peri, ~] = coe2eci(a, e, i_deg, RAAN0, om0, 0);
    scatter3(ax, r_peri(1), r_peri(2), r_peri(3), 80, 'd', ...
        'MarkerFaceColor', [1.00 0.55 0.10], 'MarkerEdgeColor', 'none', ...
        'DisplayName', sprintf('Perigee  (%.0f km alt)', orb.alt_peri));
end

%% ── Angular momentum vector (scaled arrow direction) ────────────────────────
h_hat = cross(orb.r_vec, orb.v_vec);
h_hat = h_hat / norm(h_hat) * r_max * 0.55;
quiver3(ax, 0, 0, 0, h_hat(1), h_hat(2), h_hat(3), 0, ...
    'Color', [0.70 0.40 1.00 0.80], 'LineWidth', 1.4, 'MaxHeadSize', 0.4, ...
    'DisplayName', 'h (ang. momentum)');

%% ── Labels and title ─────────────────────────────────────────────────────────
txtCol = [0.85 0.85 0.88];
xlabel(ax, 'X_{ECI} (km)', 'Color', txtCol, 'FontSize', 10);
ylabel(ax, 'Y_{ECI} (km)', 'Color', txtCol, 'FontSize', 10);
zlabel(ax, 'Z_{ECI} (km)', 'Color', txtCol, 'FontSize', 10);

title_str = sprintf('%s (ECI)  |  a = %.0f km  |  e = %.4f  |  i = %.2f°  |  T = %.1f min', ...
    upper(orb.type), a, e, orb.i, T/60);
title(ax, title_str, 'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');

legend(ax, 'Location', 'northeast', 'TextColor', txtCol, ...
    'Color', [0.08 0.10 0.15], 'EdgeColor', [0.30 0.32 0.40]);
end

%% ── Local helpers ────────────────────────────────────────────────────────────
function E = keplerSolve(M, e)
    E = M;
    for k = 1:50
        dE = (M - E + e*sin(E)) / (1 - e*cos(E));
        E  = E + dE;
        if abs(dE) < 1e-13, break; end
    end
end
