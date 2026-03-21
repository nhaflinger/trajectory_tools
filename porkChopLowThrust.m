function fig = porkChopLowThrust(departBody, arrivalBody, departJD, tofDays, options)
%PORKCHOPLOWTHRUST  Low-thrust pork-chop plot: propellant fraction vs. launch window.
%
%   fig = porkChopLowThrust(departBody, arrivalBody, departJD, tofDays)
%   fig = porkChopLowThrust(departBody, arrivalBody, departJD, tofDays, options)
%
%   For each (departure date, TOF) grid point, computes the total low-thrust
%   mission ΔV using three Edelbaum spirals (departure escape + heliocentric
%   cruise + arrival capture).  Actual planet positions are used, so orbital
%   eccentricity and phasing create the launch-window structure visible in
%   the plot.
%
%   The colour axis shows propellant mass fraction (propellant/wet mass),
%   which is independent of absolute spacecraft size and directly expresses
%   what fraction of the vehicle must be propellant.
%
%   An impulsive contour overlay (white dashed) shows the equivalent Lambert
%   ΔV from porkChopPlot for direct comparison.
%
%   Inputs:
%     departBody  - body struct from constants()
%     arrivalBody - body struct from constants()
%     departJD    - vector of departure Julian Dates
%     tofDays     - vector of TOF values (days)
%     options     (optional struct):
%       .isp                 specific impulse for propellant calc (s) [3000]
%       .wetMass             spacecraft wet mass (kg) — for prop calc  [1000]
%       .departureAltitude   km                                         [200]
%       .arrivalAltitude     km                                         [400]
%       .colorMode           'propFrac' (default) or 'deltaV'
%       .cLimPct             colorscale percentile clamp [lo hi]      [5 85]
%       .showImpulsive       overlay impulsive Lambert contours        [true]
%       .impulsiveContours   ΔV levels for impulsive overlay (km/s)   [auto]

if nargin < 5, options = struct(); end
if ~isfield(options,'isp'),               options.isp               = 3000;      end
if ~isfield(options,'wetMass'),           options.wetMass           = 1000;      end
if ~isfield(options,'departureAltitude'), options.departureAltitude = 200;       end
if ~isfield(options,'arrivalAltitude'),   options.arrivalAltitude   = 400;       end
if ~isfield(options,'colorMode'),         options.colorMode         = 'propFrac'; end
if ~isfield(options,'cLimPct'),           options.cLimPct           = [5 85];    end
if ~isfield(options,'showImpulsive'),     options.showImpulsive     = true;      end

consts = constants();
muSun  = consts.Sun.mu;
g0     = 9.80665e-3;   % km/s^2

nd = numel(departJD);
nt = numel(tofDays);

% ---- Grid computation (analytical Edelbaum only — no RK4) -----------
% Pre-compute SOI and parking-orbit ΔV constants (independent of date)
r_SOI_dep  = departBody.a  * (departBody.mu  / muSun)^(2/5);
r_SOI_arr  = arrivalBody.a * (arrivalBody.mu / muSun)^(2/5);
r_park_dep = departBody.radius  + options.departureAltitude;
r_park_arr = arrivalBody.radius + options.arrivalAltitude;

v_park_dep = sqrt(departBody.mu  / r_park_dep);
v_SOI_dep  = sqrt(departBody.mu  / r_SOI_dep);
v_SOI_arr  = sqrt(arrivalBody.mu / r_SOI_arr);
v_park_arr = sqrt(arrivalBody.mu / r_park_arr);

dv_dep_spiral = abs(v_park_dep - v_SOI_dep);   % Edelbaum dep escape
dv_arr_spiral = abs(v_SOI_arr  - v_park_arr);  % Edelbaum arr capture

% Pre-fetch actual heliocentric planet positions for all dates
r1_all = zeros(1, nd);
r2_all = zeros(nd, nt);

for i = 1:nd
    try
        r1_vec   = orbitalState(departBody, departJD(i));
        r1_all(i) = norm(r1_vec);
    catch
        r1_all(i) = departBody.a;
    end
end
for i = 1:nd
    for j = 1:nt
        try
            r2_vec      = orbitalState(arrivalBody, departJD(i) + tofDays(j));
            r2_all(i,j) = norm(r2_vec);
        catch
            r2_all(i,j) = arrivalBody.a;
        end
    end
end

% Build grids
dvGrid   = nan(nt, nd);
dvImpGrid = nan(nt, nd);

for i = 1:nd
    for j = 1:nt
        r1 = r1_all(i);
        r2 = r2_all(i, j);

        % Heliocentric Edelbaum: |v_circ(r1) - v_circ(r2)|
        v1 = sqrt(muSun / r1);
        v2 = sqrt(muSun / r2);
        dv_helio = abs(v1 - v2);

        dv_total = dv_dep_spiral + dv_helio + dv_arr_spiral;
        dvGrid(j, i) = dv_total;

        % Impulsive Lambert for overlay
        if options.showImpulsive
            try
                tof_s = tofDays(j) * 86400;
                [r1v, v1b] = orbitalState(departBody,  departJD(i));
                [r2v, v2b] = orbitalState(arrivalBody, departJD(i) + tofDays(j));
                [v1t, v2t] = lambertSolver(r1v, r2v, tof_s, muSun);
                dvImpGrid(j, i) = norm(v1t - v1b) + norm(v2b - v2t);
            catch
                dvImpGrid(j, i) = NaN;
            end
        end
    end
end

% Convert ΔV to propellant fraction via Tsiolkovsky
propFracGrid = 1 - exp(-dvGrid / (options.isp * g0));

% Select colour data
switch lower(options.colorMode)
    case 'deltav'
        cData  = dvGrid;
        cLabel = 'Total mission \DeltaV (km/s)';
    otherwise
        cData  = propFracGrid;
        cLabel = sprintf('Propellant fraction  m_{prop}/m_0   (Isp = %d s)', options.isp);
end

% Clamp colorscale
validC  = cData(isfinite(cData));
if ~isempty(validC)
    cLo = prctile(validC, options.cLimPct(1));
    cHi = prctile(validC, options.cLimPct(2));
else
    cLo = 0;  cHi = 1;
end

% ---- Figure setup ---------------------------------------------------
bgCol  = [0.06 0.06 0.10];
txtCol = [0.85 0.85 0.85];
axCol  = [0.68 0.68 0.68];

figName = sprintf('Low-Thrust Pork Chop: %s \x2192 %s', departBody.name, arrivalBody.name);
fig = figure('Name', figName, 'NumberTitle', 'off', ...
    'Color', bgCol, 'Position', [80 80 900 560]);

ax = axes('Parent', fig, ...
    'Position', [0.10 0.13 0.74 0.78], ...
    'Color', bgCol, 'XColor', axCol, 'YColor', axCol, ...
    'GridColor', [0.22 0.22 0.28], 'Box', 'on');
hold(ax, 'on');

% Filled contour
[DD, TT] = meshgrid(departJD, tofDays);
contourf(ax, DD, TT, cData, 24, 'LineColor', 'none');
colormap(ax, flipud(turbo(256)));   % turbo: low = blue (efficient), high = red (costly)
caxis(ax, [cLo cHi]);

cb = colorbar(ax, 'Location', 'eastoutside', 'Color', axCol);
cb.Label.String   = cLabel;
cb.Label.Color    = txtCol;
cb.Label.FontSize = 9;

% Impulsive Lambert contour overlay
if options.showImpulsive && any(isfinite(dvImpGrid(:)))
    validImp = dvImpGrid(isfinite(dvImpGrid));
    if isfield(options, 'impulsiveContours')
        iLevels = options.impulsiveContours;
    else
        iLevels = linspace(prctile(validImp,5), prctile(validImp,85), 6);
        iLevels = unique(round(iLevels * 10) / 10);
    end
    [C, h] = contour(ax, DD, TT, dvImpGrid, iLevels, ...
        'Color', [0.85 0.85 0.85], 'LineWidth', 0.9, 'LineStyle', '--');
    clabel(C, h, 'FontSize', 7, 'Color', [0.85 0.85 0.85], 'LabelSpacing', 300);
end

% ---- Axes labels & title --------------------------------------------
ylabel(ax, 'Time of Flight (days)', 'Color', txtCol);
title(ax, sprintf('Low-Thrust Pork Chop: %s \\rightarrow %s', ...
    departBody.name, arrivalBody.name), 'Color', txtCol, 'FontSize', 11);

% ---- Calendar ticks (identical to porkChopPlot.m) -------------------
dn_min = min(departJD) - 1721058.5;
dn_max = max(departJD) - 1721058.5;

dv_lo  = datevec(dn_min);
dv_hi  = datevec(dn_max);
yr_dns = arrayfun(@(y) datenum(y,1,1), dv_lo(1) : dv_hi(1)+1);
yr_dns = yr_dns(yr_dns >= dn_min & yr_dns <= dn_max);

y = dv_lo(1);  m = dv_lo(2) + 1;
if m > 12, m = 1; y = y + 1; end
mo_dns = [];
while datenum(y, m, 1) <= dn_max
    if m ~= 1, mo_dns(end+1) = datenum(y, m, 1); end %#ok<AGROW>
    m = m + 1;
    if m > 12, m = 1; y = y + 1; end
end

if numel(yr_dns) >= 2
    set(ax, 'XTick', yr_dns + 1721058.5, ...
            'XTickLabel', datestr(yr_dns(:), 'yyyy'));
    if ~isempty(mo_dns)
        ax.XAxis.MinorTickValues = mo_dns + 1721058.5;
        ax.XMinorTick = 'on';
    end
else
    all_dns = sort([yr_dns, mo_dns]);
    set(ax, 'XTick', all_dns + 1721058.5, ...
            'XTickLabel', datestr(all_dns(:), 'mmm yyyy'));
    xtickangle(ax, 45);
end

% Annotation: method note
noteStr = sprintf( ...
    ['3-phase Edelbaum model:\n' ...
     '  Phase 1: %s escape spiral (alt %.0f km \x2192 SOI %.0f km)\n' ...
     '  Phase 2: heliocentric spiral (%s \x2192 %s orbit)\n' ...
     '  Phase 3: %s capture spiral (SOI %.0f km \x2192 alt %.0f km)\n' ...
     'White dashed: impulsive Lambert \x0394V (km/s)\n' ...
     'Isp = %d s'], ...
    departBody.name, options.departureAltitude, r_SOI_dep, ...
    departBody.name, arrivalBody.name, ...
    arrivalBody.name, r_SOI_arr, options.arrivalAltitude, ...
    options.isp);
text(ax, 0.01, 0.01, noteStr, 'Units', 'normalized', ...
    'VerticalAlignment', 'bottom', 'Color', [txtCol 0.8], 'FontSize', 7, ...
    'BackgroundColor', [0.08 0.08 0.14 0.8], 'Margin', 4);
end
