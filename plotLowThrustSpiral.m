function fig = plotLowThrustSpiral(result, body, options)
%PLOTLOWTHRUSTSPIRALM  Visualise a low-thrust spiral orbit transfer.
%
%   fig = plotLowThrustSpiral(result, body)
%   fig = plotLowThrustSpiral(result, body, options)
%
%   Three-panel figure:
%     Left  — 2-D spiral in the orbital plane, coloured by elapsed time
%     Right-top    — altitude vs. time
%     Right-bottom — spacecraft mass vs. time
%
%   Inputs:
%     result   - struct returned by lowThrustSpiral()
%     body     - central body struct from constants()
%     options  (optional struct):
%       .title  override figure window title string
%       .units  'km' (default) or 'AU' for heliocentric spirals

if nargin < 3, options = struct(); end
if ~isfield(options, 'units'), options.units = 'km'; end

traj = result.trajectory;
if ~isfield(traj, 'x') || isempty(traj.x)
    warning('plotLowThrustSpiral: no trajectory (set thrustN > 0). Nothing to plot.');
    fig = [];
    return
end

% Unit scaling --------------------------------------------------------
switch lower(options.units)
    case 'au'
        AU   = 1.496e8;
        sc   = 1 / AU;
        uLbl = 'AU';
    otherwise
        sc   = 1;
        uLbl = 'km';
end

% Colour palette (dark theme) -----------------------------------------
bgCol  = [0.06 0.06 0.10];
txtCol = [0.85 0.85 0.85];
axCol  = [0.68 0.68 0.68];
gridC  = [0.20 0.20 0.26];

%% Figure and axes layout ---------------------------------------------
figTitle = sprintf('Low-Thrust Spiral — %s', body.name);
if isfield(options, 'title'), figTitle = options.title; end

fig = figure('Name', figTitle, 'NumberTitle', 'off', ...
    'Color', bgCol, 'Position', [80 80 1100 540]);

% Left: spiral (occupies full height)
ax1 = axes('Parent', fig, ...
    'Position', [0.055, 0.10, 0.50, 0.82], ...
    'Color', bgCol, 'XColor', axCol, 'YColor', axCol, ...
    'GridColor', gridC, 'MinorGridColor', gridC, 'Box', 'on');
grid(ax1, 'on');  hold(ax1, 'on');  axis(ax1, 'equal');

% Right-top: altitude
ax2 = axes('Parent', fig, ...
    'Position', [0.62, 0.57, 0.35, 0.35], ...
    'Color', bgCol, 'XColor', axCol, 'YColor', axCol, ...
    'GridColor', gridC, 'Box', 'on');
grid(ax2, 'on');  hold(ax2, 'on');

% Right-bottom: mass
ax3 = axes('Parent', fig, ...
    'Position', [0.62, 0.10, 0.35, 0.35], ...
    'Color', bgCol, 'XColor', axCol, 'YColor', axCol, ...
    'GridColor', gridC, 'Box', 'on');
grid(ax3, 'on');  hold(ax3, 'on');

%% Spiral plot --------------------------------------------------------
nPts  = numel(traj.t);
cmap  = parula(256);
tNorm = (traj.t - traj.t(1)) / max(traj.t(end) - traj.t(1), eps);
cidx  = max(1, min(256, floor(tNorm * 255) + 1));

for i = 1:nPts-1
    plot(ax1, traj.x([i i+1]) * sc, traj.y([i i+1]) * sc, '-', ...
        'Color', cmap(cidx(i), :), 'LineWidth', 0.9);
end

% Reference orbit circles
theta_c = linspace(0, 2*pi, 300);
r0 = result.details.r0;
r1 = result.details.r1;
plot(ax1, r0*cos(theta_c)*sc, r0*sin(theta_c)*sc, '--', ...
    'Color', [0.30 0.70 1.00 0.55], 'LineWidth', 1.0);
plot(ax1, r1*cos(theta_c)*sc, r1*sin(theta_c)*sc, '--', ...
    'Color', [1.00 0.45 0.15 0.55], 'LineWidth', 1.0);

% Start/end markers
plot(ax1, r0*sc, 0, 's', 'MarkerSize', 7, ...
    'MarkerFaceColor', [0.30 0.70 1.00], 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
plot(ax1, traj.x(end)*sc, traj.y(end)*sc, 'd', 'MarkerSize', 7, ...
    'MarkerFaceColor', [1.00 0.45 0.15], 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);

% Central body disc
bCol = bodyColor(body.name);
rb   = body.radius * sc;
th_b = linspace(0, 2*pi, 90);
fill(ax1, rb*cos(th_b), rb*sin(th_b), bCol, ...
    'EdgeColor', bCol * 0.7, 'LineWidth', 0.6, 'FaceAlpha', 0.9);

% Colorbar (time)
colormap(ax1, parula);
cb = colorbar(ax1, 'Location', 'eastoutside', 'Color', axCol);
cb.Label.String  = 'Elapsed fraction of transfer';
cb.Label.Color   = txtCol;
cb.Label.FontSize = 8;
caxis(ax1, [0 1]);

xlabel(ax1, ['x (' uLbl ')'], 'Color', txtCol);
ylabel(ax1, ['y (' uLbl ')'], 'Color', txtCol);

if ~isnan(result.tofDays)
    ttl = sprintf('%s  \\DeltaV = %.3f km/s    TOF = %.0f days', ...
        body.name, result.deltaV, result.tofDays);
else
    ttl = sprintf('%s  \\DeltaV = %.3f km/s', body.name, result.deltaV);
end
title(ax1, ttl, 'Color', txtCol, 'FontSize', 10);

% Annotation text box
if ~isnan(result.propellantMass)
    propLine = sprintf('Prop: %.0f kg  (%.1f%% wet)     m_f: %.0f kg', ...
        result.propellantMass, ...
        100 * result.propellantMass / result.details.wetMass, ...
        result.finalMass);
else
    propLine = '(no thruster — Edelbaum only)';
end
infoStr = sprintf( ...
    ['\\DeltaV (Edelbaum): %.3f km/s\n' ...
     'r_0 = %.0f km  (alt %.0f km)\n' ...
     'r_1 = %.0f km  (alt %.0f km)\n' ...
     'Isp = %.0f s    F = %.0f mN\n' ...
     '%s'], ...
    result.deltaV, ...
    r0, r0 - body.radius, ...
    r1, r1 - body.radius, ...
    result.details.isp, result.details.thrust * 1000, ...
    propLine);
text(ax1, 0.02, 0.98, infoStr, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'Color', txtCol, 'FontSize', 7.5, ...
    'FontName', 'Monospaced', ...
    'BackgroundColor', [0.10 0.10 0.16 0.85], 'Margin', 5);

%% Altitude history ---------------------------------------------------
altKm = traj.r - body.radius;
tDays = traj.t / 86400;

plot(ax2, tDays, altKm, '-', 'Color', [0.35 0.75 1.00], 'LineWidth', 1.4);
% Mark initial and target altitudes
xL = [tDays(1) tDays(end)];
plot(ax2, xL, repmat(r0 - body.radius, 1, 2), '--', ...
    'Color', [0.30 0.70 1.00 0.5], 'LineWidth', 0.8);
plot(ax2, xL, repmat(r1 - body.radius, 1, 2), '--', ...
    'Color', [1.00 0.45 0.15 0.5], 'LineWidth', 0.8);

xlabel(ax2, 'Time (days)', 'Color', txtCol, 'FontSize', 8);
ylabel(ax2, 'Altitude (km)', 'Color', txtCol, 'FontSize', 8);
title(ax2, 'Altitude Profile', 'Color', txtCol, 'FontSize', 9);
xlim(ax2, [tDays(1) tDays(end)]);

%% Mass history -------------------------------------------------------
plot(ax3, tDays, traj.mass, '-', 'Color', [1.00 0.65 0.25], 'LineWidth', 1.4);
if ~isnan(result.finalMass)
    plot(ax3, xL, repmat(result.finalMass, 1, 2), '--', ...
        'Color', [1.00 0.65 0.25 0.5], 'LineWidth', 0.8);
end

xlabel(ax3, 'Time (days)', 'Color', txtCol, 'FontSize', 8);
ylabel(ax3, 'Mass (kg)',   'Color', txtCol, 'FontSize', 8);
title(ax3, 'Spacecraft Mass', 'Color', txtCol, 'FontSize', 9);
xlim(ax3, [tDays(1) tDays(end)]);
end


%% Body colour helper (mirrors plotFlybySequence.m / tisserandGraph.m) -----
function col = bodyColor(name)
switch lower(name)
    case 'mercury',  col = [0.60 0.58 0.55];
    case 'venus',    col = [0.85 0.75 0.35];
    case 'earth',    col = [0.20 0.45 0.75];
    case 'mars',     col = [0.72 0.28 0.18];
    case 'jupiter',  col = [0.75 0.60 0.45];
    case 'saturn',   col = [0.85 0.80 0.55];
    case 'uranus',   col = [0.55 0.85 0.85];
    case 'neptune',  col = [0.25 0.40 0.85];
    case 'sun',      col = [1.00 0.85 0.20];
    case 'moon',     col = [0.75 0.75 0.75];
    otherwise,       col = [0.60 0.60 0.65];
end
end
