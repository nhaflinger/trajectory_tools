function fig = plotLowThrustInterplanetary(result, departBody, arrivalBody, options)
%PLOTLOWTHRUSTINTERPLANETARY  Visualise a three-phase low-thrust interplanetary transfer.
%
%   fig = plotLowThrustInterplanetary(result, departBody, arrivalBody)
%   fig = plotLowThrustInterplanetary(result, departBody, arrivalBody, options)
%
%   Four-panel figure:
%     Left (large) — heliocentric view: planet orbits, cruise trajectory
%                    coloured by elapsed time fraction, planet positions
%                    at departure and arrival, Sun at origin.
%     Top-right    — departure escape spiral (body-centred), or skipped note.
%     Mid-right    — arrival capture spiral  (body-centred), or skipped note.
%     Bot-right    — spacecraft mass vs. cumulative elapsed days (all phases).
%
%   Inputs:
%     result      - struct returned by lowThrustInterplanetary()
%     departBody  - body struct from constants()
%     arrivalBody - body struct from constants()
%     options     (optional struct):
%       .AU            km per AU override            [1.496e8]
%       .orbitTraceN   points per planet orbit trace [360]
%       .lvEscapeDV    ΔV (km/s) provided by launch vehicle for escape [NaN]
%                      Displayed in the departure panel when Phase 1 is skipped.
%       .lvParkAlt     parking orbit altitude (km) for LV escape annotation [200]

if nargin < 4, options = struct(); end
if ~isfield(options, 'AU'),          options.AU          = 1.496e8; end
if ~isfield(options, 'orbitTraceN'), options.orbitTraceN = 360;     end
if ~isfield(options, 'lvEscapeDV'),  options.lvEscapeDV  = NaN;     end
if ~isfield(options, 'lvParkAlt'),   options.lvParkAlt   = 200;     end

AU = options.AU;

% ---- Colour theme -------------------------------------------------------
bgCol  = [0.06 0.06 0.10];
txtCol = [0.85 0.85 0.85];
axCol  = [0.68 0.68 0.68];
gridC  = [0.20 0.20 0.26];

% ---- Extract phase results ----------------------------------------------
ph1 = result.phases.departure;
ph2 = result.phases.heliocentric;
ph3 = result.phases.arrival;

hasDep = isfield(ph1, 'trajectory') && isfield(ph1.trajectory, 'x') && ~isempty(ph1.trajectory.x);
hasArr = isfield(ph3, 'trajectory') && isfield(ph3.trajectory, 'x') && ~isempty(ph3.trajectory.x);
hasHel = isfield(ph2, 'trajectory') && isfield(ph2.trajectory, 'x') && ~isempty(ph2.trajectory.x);

det = result.details;

% ---- Figure layout ------------------------------------------------------
figName = sprintf('Low-Thrust Interplanetary: %s \x2192 %s', ...
    departBody.name, arrivalBody.name);
fig = figure('Name', figName, 'NumberTitle', 'off', ...
    'Color', bgCol, 'Position', [60 60 1250 680]);

% Main heliocentric panel (left)
axH = axes('Parent', fig, ...
    'Position', [0.04, 0.09, 0.52, 0.85], ...
    'Color', bgCol, 'XColor', axCol, 'YColor', axCol, ...
    'GridColor', gridC, 'Box', 'on');
grid(axH, 'on');  hold(axH, 'on');  axis(axH, 'equal');

% Right-top: departure spiral
axD = axes('Parent', fig, ...
    'Position', [0.60, 0.68, 0.37, 0.26], ...
    'Color', bgCol, 'XColor', axCol, 'YColor', axCol, ...
    'GridColor', gridC, 'Box', 'on');
grid(axD, 'on');  hold(axD, 'on');  axis(axD, 'equal');

% Right-mid: arrival spiral
axA = axes('Parent', fig, ...
    'Position', [0.60, 0.37, 0.37, 0.26], ...
    'Color', bgCol, 'XColor', axCol, 'YColor', axCol, ...
    'GridColor', gridC, 'Box', 'on');
grid(axA, 'on');  hold(axA, 'on');  axis(axA, 'equal');

% Right-bot: mass history
axM = axes('Parent', fig, ...
    'Position', [0.60, 0.09, 0.37, 0.22], ...
    'Color', bgCol, 'XColor', axCol, 'YColor', axCol, ...
    'GridColor', gridC, 'Box', 'on');
grid(axM, 'on');  hold(axM, 'on');

% =========================================================================
%% Heliocentric panel
% =========================================================================

theta_c = linspace(0, 2*pi, options.orbitTraceN);

% -- Planet orbit traces (circular approx at actual helio distances) -------
r1_helio = det.r1_helio;
r2_helio = det.r2_helio;

depCol = bodyColor(departBody.name);
arrCol = bodyColor(arrivalBody.name);

plot(axH, r1_helio/AU * cos(theta_c), r1_helio/AU * sin(theta_c), ...
    '-', 'Color', [depCol 0.35], 'LineWidth', 0.8);
plot(axH, r2_helio/AU * cos(theta_c), r2_helio/AU * sin(theta_c), ...
    '-', 'Color', [arrCol 0.35], 'LineWidth', 0.8);

% -- Heliocentric cruise trajectory ---------------------------------------
if hasHel
    xH = ph2.trajectory.x / AU;
    yH = ph2.trajectory.y / AU;
    tH = ph2.trajectory.t;
    nH = numel(tH);

    cmap = parula(256);
    tNorm = (tH - tH(1)) / max(tH(end) - tH(1), eps);
    cidx  = max(1, min(256, floor(tNorm * 255) + 1));
    for k = 1:nH-1
        plot(axH, xH([k k+1]), yH([k k+1]), '-', ...
            'Color', cmap(cidx(k), :), 'LineWidth', 1.2);
    end
    colormap(axH, parula);
    caxis(axH, [0 1]);
    cb = colorbar(axH, 'Location', 'southoutside', 'Color', axCol, ...
        'Position', [0.04 0.025 0.52 0.018]);
    cb.Label.String  = 'Elapsed fraction of heliocentric cruise';
    cb.Label.Color   = txtCol;
    cb.Label.FontSize = 8;
else
    % Edelbaum only — draw a straight-line arc placeholder
    ang1 = 0;
    ang2 = pi * (r2_helio > r1_helio) + pi/4;   % rough
    th_arc = linspace(ang1, ang2, 60);
    r_arc  = linspace(r1_helio, r2_helio, 60) / AU;
    plot(axH, r_arc .* cos(th_arc), r_arc .* sin(th_arc), '--', ...
        'Color', [0.70 0.70 0.90], 'LineWidth', 1.0);
    text(axH, 0, 0, 'No RK4 trajectory (thrustN = 0)', ...
        'Color', txtCol, 'FontSize', 8, 'HorizontalAlignment', 'center', ...
        'Units', 'normalized', 'Position', [0.5 0.5]);
end

% -- Planet positions at departure and arrival ----------------------------
% Departure planet at start of cruise trajectory (or angle 0 if no trajectory)
if hasHel
    depAngle = atan2(ph2.trajectory.y(1), ph2.trajectory.x(1));
else
    depAngle = 0;
end

dep_x = r1_helio/AU * cos(depAngle);
dep_y = r1_helio/AU * sin(depAngle);

% Arrival planet: last point of cruise trajectory
if hasHel
    arr_x = xH(end);
    arr_y = yH(end);
else
    arr_x = r2_helio/AU * cos(depAngle + pi * 0.7);
    arr_y = r2_helio/AU * sin(depAngle + pi * 0.7);
end

markerSz = 9;
plot(axH, dep_x, dep_y, 'o', 'MarkerSize', markerSz, ...
    'MarkerFaceColor', depCol, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
plot(axH, arr_x, arr_y, 'o', 'MarkerSize', markerSz, ...
    'MarkerFaceColor', arrCol, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);

% Labels
text(axH, dep_x, dep_y, sprintf('  %s\n  dep', departBody.name), ...
    'Color', txtCol, 'FontSize', 8, 'VerticalAlignment', 'middle');
text(axH, arr_x, arr_y, sprintf('  %s\n  arr', arrivalBody.name), ...
    'Color', txtCol, 'FontSize', 8, 'VerticalAlignment', 'middle');

% -- Sun at origin --------------------------------------------------------
sunR = max(r1_helio, r2_helio) * 0.025 / AU;   % visual radius ~2.5% of plot scale
th_s = linspace(0, 2*pi, 60);
fill(axH, sunR*cos(th_s), sunR*sin(th_s), [1.00 0.85 0.20], ...
    'EdgeColor', [1.00 0.95 0.50], 'LineWidth', 0.8, 'FaceAlpha', 0.95);

% -- Axes labels ----------------------------------------------------------
xlabel(axH, 'x (AU)', 'Color', txtCol);
ylabel(axH, 'y (AU)', 'Color', txtCol);
title(axH, sprintf('Heliocentric cruise  %s \\rightarrow %s', ...
    departBody.name, arrivalBody.name), 'Color', txtCol, 'FontSize', 10);

% -- Annotation box -------------------------------------------------------
if ~isnan(options.lvEscapeDV)
    lvLine = sprintf('  LV escape:   %.3f km/s  (not ion engine)\n', options.lvEscapeDV);
    dvTotalLine = sprintf('\\DeltaV (ion): %.3f km/s\n', result.deltaV);
else
    lvLine = '';
    dvTotalLine = sprintf('\\DeltaV total: %.3f km/s\n', result.deltaV);
end
noteStr = sprintf( ...
    ['%s' ...
     '%s' ...
     '  Departure:    %.3f km/s\n' ...
     '  Heliocentric: %.3f km/s\n' ...
     '  Arrival:      %.3f km/s\n' ...
     'Propellant:  %.1f kg  (%.1f%%)\n' ...
     'Final mass:  %.1f kg'], ...
    dvTotalLine, lvLine, ...
    det.dvDeparture, det.dvHeliocentric, det.dvArrival, ...
    result.propellantMass, ...
    100 * result.propellantMass / (result.propellantMass + result.finalMass), ...
    result.finalMass);
text(axH, 0.02, 0.98, noteStr, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'Color', txtCol, 'FontSize', 7.5, ...
    'FontName', 'Monospaced', ...
    'BackgroundColor', [0.10 0.10 0.18 0.85], 'Margin', 5);

% =========================================================================
%% Departure spiral panel
% =========================================================================
if hasDep
    plotSpiralPanel(axD, ph1, departBody, 'km', ...
        [0.30 0.70 1.00], txtCol, ...
        sprintf('%s escape spiral', departBody.name));
else
    if ~isnan(options.lvEscapeDV)
        % Draw a schematic: parking orbit circle + hyperbolic escape arc
        r_park = departBody.radius + options.lvParkAlt;
        % Parking orbit circle
        th_c = linspace(0, 2*pi, 200);
        plot(axD, r_park*cos(th_c), r_park*sin(th_c), '-', ...
            'Color', [0.30 0.70 1.00 0.5], 'LineWidth', 1.0);
        % Hyperbolic escape arc: r = p/(1 + e*cos(theta)), e=2, periapsis at r_park
        % p = r_park*(e^2-1) = r_park*3; at theta=0: r = p/(1+e) = r_park ✓
        e_hyp  = 2.0;
        p_hyp  = r_park * (e_hyp^2 - 1);
        th_esc = linspace(-1.1, 1.1, 100);   % asymptote half-angle ~acos(-1/e)=2.09 rad
        r_esc  = p_hyp ./ (1 + e_hyp * cos(th_esc));
        rot    = pi/4;   % rotate 45° so arc opens upper-right
        x_esc  = r_esc .* cos(th_esc + rot);
        y_esc  = r_esc .* sin(th_esc + rot);
        plot(axD, x_esc, y_esc, '-', 'Color', [1.00 0.82 0.20], 'LineWidth', 1.5);
        plot(axD, x_esc(end), y_esc(end), '>', 'MarkerSize', 7, ...
            'MarkerFaceColor', [1.00 0.82 0.20], 'MarkerEdgeColor', 'w');
        % Central body disc
        bCol = bodyColor(departBody.name);
        rb   = departBody.radius;
        fill(axD, rb*cos(th_c), rb*sin(th_c), bCol, ...
            'EdgeColor', bCol*0.7, 'LineWidth', 0.6, 'FaceAlpha', 0.9);
        axis(axD, 'equal');
        % ΔV annotation
        dvStr = sprintf( ...
            ['LV provides escape\n' ...
             '\\DeltaV = %.3f km/s\n' ...
             'from %.0f km alt\n' ...
             '(ion engine off)'], ...
            options.lvEscapeDV, options.lvParkAlt);
        text(axD, 0.05, 0.97, dvStr, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'Color', [1.00 0.90 0.50], 'FontSize', 8, ...
            'FontWeight', 'bold', ...
            'BackgroundColor', [0.10 0.10 0.16 0.85], 'Margin', 4);
    else
        text(axD, 0.5, 0.5, sprintf('Phase 1 skipped\n(LV provides escape)'), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', [txtCol 0.7], 'FontSize', 9);
    end
    title(axD, sprintf('%s escape — launch vehicle', departBody.name), ...
        'Color', txtCol, 'FontSize', 9);
end

% =========================================================================
%% Arrival spiral panel
% =========================================================================
if hasArr
    plotSpiralPanel(axA, ph3, arrivalBody, 'km', ...
        [1.00 0.45 0.15], txtCol, ...
        sprintf('%s capture spiral', arrivalBody.name));
else
    axis(axA, 'off');
    text(axA, 0.5, 0.5, sprintf('Phase 3 skipped\n(no arrival capture)'), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', [txtCol 0.7], 'FontSize', 9);
    title(axA, sprintf('%s capture spiral', arrivalBody.name), ...
        'Color', txtCol, 'FontSize', 9);
end

% =========================================================================
%% Mass history (all phases concatenated)
% =========================================================================
tOffset = 0;
phaseColors = {[0.30 0.70 1.00], [0.55 0.90 0.55], [1.00 0.45 0.15]};
phaseNames  = {'Departure escape', 'Heliocentric cruise', 'Arrival capture'};
phases      = {ph1, ph2, ph3};
hasPhase    = {hasDep, hasHel, hasArr};

legendH = [];
legendL = {};

for k = 1:3
    ph = phases{k};
    if hasPhase{k}
        tDays = ph.trajectory.t / 86400 + tOffset;
        mass  = ph.trajectory.mass;
        h = plot(axM, tDays, mass, '-', 'Color', phaseColors{k}, 'LineWidth', 1.5);
        legendH(end+1) = h; %#ok<AGROW>
        legendL{end+1} = phaseNames{k}; %#ok<AGROW>
        tOffset = tOffset + ph.trajectory.t(end) / 86400;
    elseif ~isnan(ph.tofDays) && ph.tofDays > 0
        % Has analytical TOF but no trajectory — draw a line between masses
        % (approximate mass drawdown)
    end
end

if ~isempty(legendH)
    lg = legend(axM, legendH, legendL, 'Location', 'northeast', ...
        'TextColor', txtCol, 'FontSize', 7);
    lg.Color = [0.10 0.10 0.16];
    lg.EdgeColor = axCol;
end

xlabel(axM, 'Elapsed mission time (days)', 'Color', txtCol, 'FontSize', 8);
ylabel(axM, 'Mass (kg)', 'Color', txtCol, 'FontSize', 8);
title(axM, 'Spacecraft mass (all phases)', 'Color', txtCol, 'FontSize', 9);
end


% =========================================================================
%% Local helper: draw a compact spiral panel
% =========================================================================
function plotSpiralPanel(ax, ph, body, units, spiralCol, txtCol, titleStr)

traj = ph.trajectory;

switch lower(units)
    case 'au',  sc = 1 / 1.496e8;  uLbl = 'AU';
    otherwise,  sc = 1;             uLbl = 'km';
end

% Spiral trajectory coloured by time
nPts  = numel(traj.t);
cmap  = parula(256);
tNorm = (traj.t - traj.t(1)) / max(traj.t(end) - traj.t(1), eps);
cidx  = max(1, min(256, floor(tNorm * 255) + 1));
for i = 1:nPts-1
    plot(ax, traj.x([i i+1])*sc, traj.y([i i+1])*sc, '-', ...
        'Color', cmap(cidx(i), :), 'LineWidth', 0.8);
end

% Initial / target orbit circles
theta_c = linspace(0, 2*pi, 200);
r0 = ph.details.r0;
r1 = ph.details.r1;
plot(ax, r0*cos(theta_c)*sc, r0*sin(theta_c)*sc, '--', ...
    'Color', [0.30 0.70 1.00 0.50], 'LineWidth', 0.8);
plot(ax, r1*cos(theta_c)*sc, r1*sin(theta_c)*sc, '--', ...
    'Color', [1.00 0.45 0.15 0.50], 'LineWidth', 0.8);

% Start / end markers
plot(ax, r0*sc, 0, 's', 'MarkerSize', 6, ...
    'MarkerFaceColor', [0.30 0.70 1.00], 'MarkerEdgeColor', 'w');
plot(ax, traj.x(end)*sc, traj.y(end)*sc, 'd', 'MarkerSize', 6, ...
    'MarkerFaceColor', [1.00 0.45 0.15], 'MarkerEdgeColor', 'w');

% Central body
bCol = bodyColor(body.name);
rb   = body.radius * sc;
th_b = linspace(0, 2*pi, 60);
fill(ax, rb*cos(th_b), rb*sin(th_b), bCol, ...
    'EdgeColor', bCol*0.7, 'LineWidth', 0.6, 'FaceAlpha', 0.9);

% Annotation
if ~isnan(ph.propellantMass)
    dvStr = sprintf('\\DeltaV = %.3f km/s\nTOF = %.0f d\nProp = %.1f kg', ...
        ph.deltaV, ph.tofDays, ph.propellantMass);
else
    dvStr = sprintf('\\DeltaV = %.3f km/s', ph.deltaV);
end
text(ax, 0.03, 0.97, dvStr, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'Color', txtCol, 'FontSize', 7, ...
    'BackgroundColor', [0.10 0.10 0.16 0.85], 'Margin', 3);

xlabel(ax, ['x (' uLbl ')'], 'Color', txtCol, 'FontSize', 7);
ylabel(ax, ['y (' uLbl ')'], 'Color', txtCol, 'FontSize', 7);
title(ax, titleStr, 'Color', txtCol, 'FontSize', 9);
end


%% Body colour helper
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
    case 'vesta',    col = [0.65 0.60 0.55];
    case 'ceres',    col = [0.55 0.55 0.60];
    otherwise,       col = [0.60 0.60 0.65];
end
end
