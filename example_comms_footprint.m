%EXAMPLE_COMMS_FOOTPRINT  Demonstration of linkBudget and sensorFootprint tools.
%
%   Sections:
%     1. Static X-band downlink budget
%     2. Static UHF command uplink budget
%     3. Pass analysis — ISS-like orbit over Austin TX
%     4. Sensor footprint — nadir imager at three altitudes
%     5. Sensor footprint — off-nadir pointing trade (ISS orbit)
%     6. Coverage swath width table

clear;  clc;

%% ── Output directory ────────────────────────────────────────────────────────
outDir = fullfile(fileparts(mfilename('fullpath')), 'output_comms_footprint');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];

fprintf('============================================================\n');
fprintf('  COMMUNICATIONS & FOOTPRINT ANALYSIS EXAMPLES\n');
fprintf('============================================================\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% Section 1: Static link budget — X-band downlink
fprintf('--- Section 1: X-band downlink (static) ---\n');

[result_xband, fig1] = linkBudget( ...
    'Freq_GHz',     8.0,    ...
    'P_tx_dBW',     3.0,    ...
    'G_tx_dBi',     6.0,    ...
    'G_rx_dBi',     48.0,   ...
    'T_sys_K',      100,    ...
    'DataRate_bps', 10e6,   ...
    'ReqEbN0_dB',   10.0,   ...
    'Range_km',     1500,   ...
    'Plot',         true);

saveas(fig1, fullfile(outDir, 'linkbudget_xband_static.png'));
fprintf('Saved: linkbudget_xband_static.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% Section 2: Static link budget — UHF command uplink
fprintf('--- Section 2: UHF command uplink (static) ---\n');

[result_uhf, fig2] = linkBudget( ...
    'Freq_GHz',     0.437,  ...
    'P_tx_dBW',     10.0,   ...
    'G_tx_dBi',     6.0,    ...
    'G_rx_dBi',     0.0,    ...
    'T_sys_K',      500,    ...
    'DataRate_bps', 9600,   ...
    'ReqEbN0_dB',   10.0,   ...
    'Range_km',     800,    ...
    'Plot',         true);

saveas(fig2, fullfile(outDir, 'linkbudget_uhf_static.png'));
fprintf('Saved: linkbudget_uhf_static.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% Section 3: Pass analysis — ISS-like orbit over Austin TX
fprintf('--- Section 3: ISS pass analysis over Austin TX ---\n');

orb_iss   = earthOrbit('circular', 410, 51.6);
gs_austin = struct('lat', 30.27, 'lon', -97.74, 'alt', 0.15, 'name', 'Austin TX');

[result_pass, fig3] = linkBudget(orb_iss, gs_austin, ...
    'Freq_GHz',     2.2,    ...
    'P_tx_dBW',     0.0,    ...
    'G_tx_dBi',     3.0,    ...
    'G_rx_dBi',     35.0,   ...
    'T_sys_K',      200,    ...
    'DataRate_bps', 1e6,    ...
    'ReqEbN0_dB',   10.0,   ...
    'Plot',         true);

saveas(fig3, fullfile(outDir, 'linkbudget_pass_ISS_Austin.png'));
fprintf('Saved: linkbudget_pass_ISS_Austin.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% Section 4: Sensor footprint — nadir-pointing imager, three altitudes
fprintf('--- Section 4: Nadir footprint vs altitude (30 deg half-angle) ---\n');

alts_km  = [400, 600, 1000];
cols4    = {[0.30 0.75 0.93], [0.25 0.90 0.55], [1.00 0.65 0.25]};
half_ang = 30;

% Create figure and axes with coastlines drawn once
fig4   = figure('Color', bgCol, 'Position', [100 80 1100 580]);
ax4    = axes('Parent', fig4, 'Color', [0.14 0.18 0.26], ...
    'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', [0.25 0.28 0.35], 'GridAlpha', 0.5);
hold(ax4, 'on');  grid(ax4, 'on');

% Coastlines
[clon4, clat4] = loadCoastlines_local();
if ~isempty(clon4)
    plot(ax4, clon4, clat4, '-', 'Color', [0.42 0.52 0.60], 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
end

for ki = 1:numel(alts_km)
    alt_i = alts_km(ki);
    col_i = cols4{ki};

    fp_i = sensorFootprint(0, 0, alt_i, half_ang, ...
        'PlotAxes',   ax4,    ...
        'FillColor',  col_i,  ...
        'FillAlpha',  0.20,   ...
        'ShowSubSat', ki == 1, ...
        'ShowCenter', false);

    % Add a legend proxy
    fill(ax4, NaN, NaN, col_i, 'FaceAlpha', 0.4, 'EdgeColor', col_i, ...
        'DisplayName', sprintf('%d km alt  (r=%.0f km, MinEl=%.1f°)', ...
        alt_i, fp_i.footprint_radius_km, fp_i.min_elevation_deg));

    fprintf('  Alt %d km: footprint radius = %.0f km, min elevation = %.1f deg\n', ...
        alt_i, fp_i.footprint_radius_km, fp_i.min_elevation_deg);
end

% Sub-satellite point
scatter(ax4, 0, 0, 100, 'x', 'MarkerEdgeColor', [0.95 0.90 0.25], 'LineWidth', 2.0, ...
    'DisplayName', 'Sub-satellite (0°N, 0°E)');

xlim(ax4, [-180 180]);  ylim(ax4, [-90 90]);
xticks(ax4, -180:30:180);  yticks(ax4, -90:30:90);
xlabel(ax4, 'Longitude (deg)', 'Color', txtCol, 'FontSize', 11);
ylabel(ax4, 'Latitude (deg)',  'Color', txtCol, 'FontSize', 11);
title(ax4, sprintf('Nadir Footprint vs Altitude  |  Half-angle = %d°  |  Sub-sat: 0°N, 0°E', half_ang), ...
    'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');
legend(ax4, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'TextColor', txtCol, 'Color', [0.08 0.10 0.14], 'EdgeColor', [0.30 0.32 0.38]);

saveas(fig4, fullfile(outDir, 'footprint_altitude_comparison.png'));
fprintf('Saved: footprint_altitude_comparison.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% Section 5: Sensor footprint — off-nadir pointing trade
fprintf('--- Section 5: Off-nadir pointing trade (ISS 410 km, 15 deg half-angle) ---\n');

point_els  = [0, 15, 30];
cols5      = {[0.30 0.75 0.93], [0.25 0.90 0.55], [1.00 0.55 0.30]};
half_ang5  = 15;
lat5       = 40;   lon5 = 0;
alt5       = 410;  pointAz5 = 90;   % pointing East

fig5   = figure('Color', bgCol, 'Position', [100 80 1100 580]);
ax5    = axes('Parent', fig5, 'Color', [0.14 0.18 0.26], ...
    'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', [0.25 0.28 0.35], 'GridAlpha', 0.5);
hold(ax5, 'on');  grid(ax5, 'on');

[clon5, clat5] = loadCoastlines_local();
if ~isempty(clon5)
    plot(ax5, clon5, clat5, '-', 'Color', [0.42 0.52 0.60], 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
end

for pi = 1:numel(point_els)
    pel   = point_els(pi);
    col_i = cols5{pi};

    fp_i = sensorFootprint(lat5, lon5, alt5, half_ang5, ...
        'PointEl_deg', pel,        ...
        'PointAz_deg', pointAz5,   ...
        'PlotAxes',    ax5,        ...
        'FillColor',   col_i,      ...
        'FillAlpha',   0.25,       ...
        'ShowSubSat',  pi == 1,    ...
        'ShowCenter',  pel > 0);

    lbl = sprintf('PointEl=%d°  (MinEl=%.1f°)', pel, fp_i.min_elevation_deg);
    fill(ax5, NaN, NaN, col_i, 'FaceAlpha', 0.4, 'EdgeColor', col_i, ...
        'DisplayName', lbl);

    fprintf('  PointEl = %d deg: footprint center = (%.2f N, %.2f E), min elevation = %.1f deg\n', ...
        pel, fp_i.lat_center, fp_i.lon_center, fp_i.min_elevation_deg);
end

scatter(ax5, lon5, lat5, 100, 'x', 'MarkerEdgeColor', [0.95 0.90 0.25], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Sub-sat (%d°N, %d°E)', lat5, lon5));

xlim(ax5, [-60 60]);   ylim(ax5, [10 70]);
xlabel(ax5, 'Longitude (deg)', 'Color', txtCol, 'FontSize', 11);
ylabel(ax5, 'Latitude (deg)',  'Color', txtCol, 'FontSize', 11);
title(ax5, sprintf('Off-Nadir Pointing Trade  |  Alt = %d km  |  Half-angle = %d°  |  Az = %d° (East)', ...
    alt5, half_ang5, pointAz5), 'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');
legend(ax5, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'TextColor', txtCol, 'Color', [0.08 0.10 0.14], 'EdgeColor', [0.30 0.32 0.38]);

saveas(fig5, fullfile(outDir, 'footprint_offnadir_trade.png'));
fprintf('Saved: footprint_offnadir_trade.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% Section 6: Coverage swath width table
fprintf('--- Section 6: Footprint radius & minimum elevation table ---\n\n');

alt_list = [400, 500, 600, 700, 800];
ha_list  = [10, 20, 30, 45];

% Header
fprintf('  Half-angle: ');
for hi = 1:numel(ha_list)
    fprintf('%3d deg         ', ha_list(hi));
end
fprintf('\n');

fprintf('  Alt (km)  |');
for hi = 1:numel(ha_list)
    fprintf(' R(km)  MinEl |');
end
fprintf('\n');
fprintf('  %s\n', repmat('-', 1, 10 + numel(ha_list)*16));

for ai = 1:numel(alt_list)
    alt_i = alt_list(ai);
    fprintf('  %-8d  |', alt_i);
    for hi = 1:numel(ha_list)
        fp_t = sensorFootprint(0, 0, alt_i, ha_list(hi));
        fprintf(' %5.0f %6.1f° |', fp_t.footprint_radius_km, fp_t.min_elevation_deg);
    end
    fprintf('\n');
end
fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════════
%% Section 7: Access Window Gantt — ISS over global ground stations
fprintf('\n--- Section 7: Access Window Gantt ---\n');

orb_gantt = earthOrbit('circular', 410, 51.6);

gs_network = struct(...
    'name',  {'Svalbard', 'Fairbanks', 'Austin TX', 'Madrid', 'Singapore', 'Santiago', 'Canberra'}, ...
    'lat',   { 78.2,       64.8,        30.3,        40.4,     1.3,        -33.4,      -35.3}, ...
    'lon',   { 15.4,      -147.7,      -97.7,        -3.7,   103.8,        -70.6,      149.1}, ...
    'alt',   {  0.5,        0.1,         0.15,         0.7,    0.015,         0.5,       0.8});

[fig7, summary7] = plotAccessWindowGantt(orb_gantt, gs_network, 24, ...
    'MinElevation', 10, ...
    'StepSize',     30, ...
    'ColorByElevation', true, ...
    'ShowMaxEl',    true);

saveas(fig7, fullfile(outDir, 'access_window_gantt.png'));
fprintf('  Saved: access_window_gantt.png\n');

% Print total contact per station
fprintf('\n  24-hour contact summary:\n');
for k = 1:numel(summary7)
    fprintf('    %-12s  %d passes  %.1f min total\n', ...
        summary7(k).name, summary7(k).n_passes, summary7(k).total_contact_min);
end

fprintf('\n============================================================\n');
fprintf('  All outputs saved to: %s\n', outDir);
fprintf('============================================================\n');

%% ── Local helper: coastline loader ─────────────────────────────────────────
function [lon, lat] = loadCoastlines_local()
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
    try
        downloadCoastlines();
        s = load(cachePath, 'coastlon', 'coastlat');
        lon = s.coastlon;  lat = s.coastlat;
    catch ME
        warning('example_comms_footprint:noCoastlines', ...
            'Could not load coastlines: %s', ME.message);
    end
end
