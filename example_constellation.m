%EXAMPLE_CONSTELLATION  Demonstrate Walker constellation design tools.
%
%   Sections:
%     1. Starlink-like shell  (24/6/1 at 550 km, 53 deg)
%     2. GPS-like MEO         (24/6/2 at 20200 km, 55 deg)
%     3. Coverage comparison  (1-plane vs 4-plane vs 6-plane at 550 km)
%     4. Repeating ground track survey  (i=98 deg)
%     5. Repeating ground track plot   (closest to 550 km)

clear; close all; clc;

%% ── Output directory ────────────────────────────────────────────────────────
outDir = fullfile(fileparts(mfilename('fullpath')), 'output_constellation');
if ~isfolder(outDir)
    mkdir(outDir);
end

bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];

fprintf('===================================================\n');
fprintf('  Walker Constellation Design Tools — Demo\n');
fprintf('===================================================\n\n');

%% ══════════════════════════════════════════════════════════════════════════════
%% Section 1: Starlink-like shell  (24/6/1 at 550 km, 53 deg)
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('--- Section 1: Starlink-like shell (Walker 24/6/1) ---\n');

sats_sl  = walkerConstellation(53, 550, 24, 6, 1);
fig1     = plotConstellationGroundTrack(sats_sl, 'NumOrbits', 1, ...
    'Title', 'Starlink-like Shell  |  Walker 24/6/1  |  550 km  |  i = 53°  |  24 satellites');

fig1_path = fullfile(outDir, 'starlink_like_groundtrack.png');
exportgraphics(fig1, fig1_path, 'Resolution', 150, 'BackgroundColor', bgCol);
fprintf('  Saved: %s\n\n', fig1_path);

%% ══════════════════════════════════════════════════════════════════════════════
%% Section 2: GPS-like MEO  (24/6/2 at 20200 km, 55 deg)
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('--- Section 2: GPS-like MEO (Walker 24/6/2) ---\n');

sats_gps = walkerConstellation(55, 20200, 24, 6, 2);
fig2     = plotConstellationGroundTrack(sats_gps, 'NumOrbits', 1, ...
    'Title', 'GPS-like MEO  |  Walker 24/6/2  |  20200 km  |  i = 55°  |  24 satellites');

fig2_path = fullfile(outDir, 'gps_like_groundtrack.png');
exportgraphics(fig2, fig2_path, 'Resolution', 150, 'BackgroundColor', bgCol);
fprintf('  Saved: %s\n\n', fig2_path);

%% ══════════════════════════════════════════════════════════════════════════════
%% Section 3: Coverage comparison
%%   Three constellations at 550 km, 53 deg, 24 sats — varying number of planes
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('--- Section 3: Coverage comparison ---\n');
fprintf('  Configurations: 24/1/0, 24/4/1, 24/6/1 at 550 km, 53 deg\n');
fprintf('  Duration: 24 hr, GridRes: 5 deg, MinElev: 10 deg, StepSize: 120 s\n\n');

configs = {
    24, 1, 0, '24/1/0 (1 plane)';
    24, 4, 1, '24/4/1 (4 planes)';
    24, 6, 1, '24/6/1 (6 planes)';
};

cov_results = cell(size(configs, 1), 1);

for c = 1:size(configs, 1)
    T_c = configs{c,1};
    P_c = configs{c,2};
    F_c = configs{c,3};
    lbl = configs{c,4};

    sats_c = walkerConstellation(53, 550, T_c, P_c, F_c);
    fprintf('  Running coverage for %s ...\n', lbl);
    cov_c  = constellationCoverage(sats_c, 24, ...
        'GridRes', 5, 'MinElevation', 10, 'StepSize', 120);
    cov_results{c} = cov_c;
end

% Print comparison table
fprintf('\n  %-20s  %8s  %10s  %12s\n', 'Config', 'Mean Cov', 'Max RevisitHr', 'Mean RevisitHr');
fprintf('  %s\n', repmat('-', 1, 60));
for c = 1:size(configs, 1)
    lbl   = configs{c,4};
    cov_c = cov_results{c};
    mean_cov   = mean(cov_c.coverage_frac(:), 'omitnan');
    max_rev    = max(cov_c.revisit_max_hr(:), [], 'omitnan');
    mean_rev   = mean(cov_c.revisit_mean_hr(~isnan(cov_c.revisit_mean_hr(:))), 'omitnan');
    fprintf('  %-20s  %8.3f  %10.2f  %12.2f\n', lbl, mean_cov, max_rev, mean_rev);
end
fprintf('\n');

% Save coverage map for 6-plane case
fig3 = plotCoverage(cov_results{3});
title(findobj(fig3, 'Type', 'axes', '-not', 'Tag', 'legend'), ...
    'Coverage: Walker 24/6/1  |  550 km  |  i = 53°  |  24 satellites  |  MinEl = 10°', ...
    'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');
fig3_path = fullfile(outDir, 'coverage_6plane_550km.png');
exportgraphics(fig3, fig3_path, 'Resolution', 150, 'BackgroundColor', bgCol);
fprintf('  Coverage map saved: %s\n\n', fig3_path);

%% ══════════════════════════════════════════════════════════════════════════════
%% Section 4: Repeating ground track survey at i=98 deg
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('--- Section 4: Repeating Ground Track Survey (i=98 deg) ---\n\n');

rgt_configs = [
    14, 1;    % 14 rev/day
    15, 1;    % 15 rev/day
    29, 2;    % 14.5 rev/day
    43, 3;    % ~14.33 rev/day
    59, 4;    % ~14.75 rev/day
];

inc_rgt = 98;   % deg

rgt_results = struct([]);
for r = 1:size(rgt_configs, 1)
    N_rev_r = rgt_configs(r,1);
    N_day_r = rgt_configs(r,2);
    res_r   = repeatingGroundTrack(N_rev_r, N_day_r, inc_rgt);
    rgt_results(r).res = res_r;
end

% Print table
fprintf('  %8s  %10s  %10s  %16s  %14s\n', ...
    'N_rev/N_day', 'Alt (km)', 'Period (min)', 'RAAN drift (deg/day)', 'Spacing (deg)');
fprintf('  %s\n', repmat('-', 1, 70));
for r = 1:size(rgt_configs, 1)
    res_r   = rgt_results(r).res;
    ratio_s = sprintf('%d/%d', res_r.N_rev, res_r.N_day);
    fprintf('  %8s  %10.2f  %10.4f  %16.6f  %14.4f\n', ...
        ratio_s, res_r.alt_km, res_r.period_s/60, ...
        res_r.RAAN_dot_deg_day, res_r.ground_track_spacing_deg);
end
fprintf('\n');

%% ══════════════════════════════════════════════════════════════════════════════
%% Section 5: Repeating ground track plot
%%   Pick the configuration closest to 550 km from Section 4
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('--- Section 5: Repeating Ground Track Plot ---\n');

% Find closest to 550 km
alts = arrayfun(@(x) x.res.alt_km, rgt_results);
[~, best_idx] = min(abs(alts - 550));
best_res = rgt_results(best_idx).res;

fprintf('  Closest to 550 km: %d/%d, alt=%.2f km\n', ...
    best_res.N_rev, best_res.N_day, best_res.alt_km);

% Create earthOrbit struct at this altitude/inclination
orb_rgt = earthOrbit('coe', best_res.a_km, best_res.e, best_res.i_deg, 0, 0, 0);

% Plot for exactly N_day full days to show repeat closure
duration_orbits = best_res.N_rev;   % N_rev orbits = N_day sidereal days exactly

fig5 = plotGroundTrack(orb_rgt, ...
    'NumOrbits',    duration_orbits, ...
    'J2',           true, ...
    'ColorByOrbit', false, ...
    'ShowNodes',    false);

% Update title to show repeating nature
ax5 = fig5.CurrentAxes;
title(ax5, sprintf('Repeating Ground Track %d/%d | alt=%.0f km | i=%.1f deg | T=%.4f min', ...
    best_res.N_rev, best_res.N_day, best_res.alt_km, best_res.i_deg, best_res.period_s/60), ...
    'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');

fig5_path = fullfile(outDir, 'repeating_groundtrack.png');
exportgraphics(fig5, fig5_path, 'Resolution', 150, 'BackgroundColor', bgCol);
fprintf('  Saved: %s\n\n', fig5_path);

fprintf('===================================================\n');
fprintf('  example_constellation.m complete.\n');
fprintf('  Outputs in: %s\n', outDir);
fprintf('===================================================\n');
