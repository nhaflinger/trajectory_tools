%% example_earth_coverage.m
% Demonstrates propagateOrbit, accessWindows, coverageAnalysis, and plotCoverage.
%
% Functions used: earthOrbit, propagateOrbit, accessWindows, coverageAnalysis,
%                 plotCoverage, plotGroundTrack

SCRIPT_DIR = fileparts(mfilename('fullpath'));
OUT_DIR    = fullfile(SCRIPT_DIR, 'output_earth_coverage');
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

fprintf('=======================================================\n');
fprintf(' Earth Coverage Analysis Example\n');
fprintf('=======================================================\n\n');

%% ── 1. Propagate ISS-like LEO for 3 days using J2 method ───────────────────
fprintf('--- 1. ISS-like orbit propagation (J2, 3 days) ---\n');

iss = earthOrbit('circular', 410, 51.6);
fprintf('  Orbit:    alt = %.0f km, i = %.1f deg, T = %.2f min\n', ...
    iss.alt_peri, iss.i, iss.period/60);

dur_3day = 3 * 86400;   % 3 days in seconds

traj_j2 = propagateOrbit(iss, dur_3day, 'Method', 'j2', 'StepSize', 30);

fprintf('  Steps:    %d (dt = 30 s)\n', numel(traj_j2.t));
fprintf('  Initial:  alt = %.2f km, lat = %.2f deg, lon = %.2f deg\n', ...
    traj_j2.alt(1), traj_j2.lat(1), traj_j2.lon(1));
fprintf('  Final:    alt = %.2f km, lat = %.2f deg, lon = %.2f deg\n', ...
    traj_j2.alt(end), traj_j2.lat(end), traj_j2.lon(end));
fprintf('  Mean alt: %.3f km,  Std: %.4f km\n', ...
    mean(traj_j2.alt), std(traj_j2.alt));
fprintf('  Alt range: %.3f to %.3f km\n\n', ...
    min(traj_j2.alt), max(traj_j2.alt));

% Ground track for 3 orbits
fig_gt = plotGroundTrack(iss, 'NumOrbits', 3, 'J2', true);
saveas(fig_gt, fullfile(OUT_DIR, 'iss_groundtrack.png'));
fprintf('  Saved: iss_groundtrack.png\n\n');

%% ── 2. Drag propagation — compare altitude decay ────────────────────────────
fprintf('--- 2. ISS orbit with drag (CdAm = 0.01 m^2/kg, 1 day) ---\n');

dur_1day = 86400;   % 1 day

traj_drag = propagateOrbit(iss, dur_1day, 'Method', 'drag', ...
    'StepSize', 60, 'CdAm', 0.01);
traj_j2_1d = propagateOrbit(iss, dur_1day, 'Method', 'j2', 'StepSize', 60);

alt_decay = traj_j2_1d.alt(end) - traj_drag.alt(end);

fprintf('  J2 only  — final alt: %.4f km\n', traj_j2_1d.alt(end));
fprintf('  J2+drag  — final alt: %.4f km\n', traj_drag.alt(end));
fprintf('  Altitude decay due to drag over 1 day: %.4f km\n', alt_decay);
fprintf('  Mean drag alt: %.4f km (vs J2 start: %.4f km)\n\n', ...
    mean(traj_drag.alt), traj_drag.alt(1));

%% ── 3. Access windows for ISS over Goddard Space Flight Center ──────────────
fprintf('--- 3. Access windows: ISS over Goddard (lat=38, lon=-77) ---\n');

goddard_lat = 38.0;
goddard_lon = -77.0;
dur_3day    = 3 * 86400;

wins = accessWindows(iss, goddard_lat, goddard_lon, dur_3day, ...
    'MinElev', 5, 'Method', 'j2');

if isempty(wins)
    fprintf('  No passes found over Goddard in %.0f hours.\n\n', dur_3day/3600);
else
    fprintf('  Passes over Goddard in %.0f hours: %d\n', dur_3day/3600, numel(wins));
    total_access = sum([wins.duration_s]);
    fprintf('  Total access time:   %.1f s  (%.2f min)\n', ...
        total_access, total_access/60);
    fprintf('  Mean pass duration:  %.1f s\n', total_access / numel(wins));
    fprintf('  Max elevation (all passes): %.2f deg\n', max([wins.max_elev_deg]));
    fprintf('\n  Pass listing (first 5):\n');
    for w = 1:min(5, numel(wins))
        fprintf('    Pass %2d:  start = %8.1f s, dur = %6.1f s, max_el = %5.2f deg\n', ...
            w, wins(w).start_s, wins(w).duration_s, wins(w).max_elev_deg);
    end
end
fprintf('\n');

%% ── 4. Coverage analysis: SSO at 550 km over 1 day ─────────────────────────
fprintf('--- 4. Coverage analysis: SSO 550 km, 1 day, GridRes=5 deg ---\n');

sso = earthOrbit('sso', 550);
fprintf('  SSO:  i = %.4f deg, alt = %.0f km, T = %.2f min\n', ...
    sso.i, sso.alt_peri, sso.period/60);

cov_sso = coverageAnalysis(sso, dur_1day, ...
    'GridRes',  5,    ...
    'MinElev',  5,    ...
    'Method',  'j2', ...
    'StepSize', 30);

valid_cov  = cov_sso.coverage_frac(cov_sso.coverage_frac > 0);
mean_cov   = mean(cov_sso.coverage_frac(:));
pct_any    = 100 * sum(cov_sso.coverage_frac(:) > 0) / numel(cov_sso.coverage_frac);
revisit_ok = cov_sso.revisit_mean_hr(~isnan(cov_sso.revisit_mean_hr));

fprintf('  Grid points:         %d x %d = %d total\n', ...
    numel(cov_sso.lat_vec), numel(cov_sso.lon_vec), ...
    numel(cov_sso.lat_vec)*numel(cov_sso.lon_vec));
fprintf('  Mean coverage frac:  %.4f\n', mean_cov);
fprintf('  Fraction with any coverage: %.1f %%\n', pct_any);
if ~isempty(revisit_ok)
    fprintf('  Mean revisit time (covered pts): %.2f hr\n', mean(revisit_ok));
    fprintf('  Max revisit time (covered pts):  %.2f hr\n', max(revisit_ok));
end
fprintf('  Max passes at any point:   %d\n\n', max(cov_sso.n_passes(:)));

fig_sso_cov = plotCoverage(cov_sso, 'Quantity', 'coverage');
saveas(fig_sso_cov, fullfile(OUT_DIR, 'sso_coverage_frac.png'));

fig_sso_rev = plotCoverage(cov_sso, 'Quantity', 'revisit');
saveas(fig_sso_rev, fullfile(OUT_DIR, 'sso_revisit_time.png'));

fprintf('  Saved: sso_coverage_frac.png, sso_revisit_time.png\n\n');

%% ── 5. Coverage analysis: ISS orbit over 1 day ──────────────────────────────
fprintf('--- 5. Coverage analysis: ISS orbit (51.6 deg), 1 day, GridRes=5 deg ---\n');

% Ground station to mark on the plot
goddard_gs = struct('lat', goddard_lat, 'lon', goddard_lon, 'name', 'Goddard');

cov_iss = coverageAnalysis(iss, dur_1day, ...
    'GridRes',  5,    ...
    'MinElev',  5,    ...
    'Method',  'j2', ...
    'StepSize', 30);

mean_cov_iss  = mean(cov_iss.coverage_frac(:));
pct_any_iss   = 100 * sum(cov_iss.coverage_frac(:) > 0) / numel(cov_iss.coverage_frac);
revisit_iss   = cov_iss.revisit_mean_hr(~isnan(cov_iss.revisit_mean_hr));

fprintf('  Mean coverage frac:  %.4f\n', mean_cov_iss);
fprintf('  Fraction with any coverage: %.1f %%\n', pct_any_iss);
if ~isempty(revisit_iss)
    fprintf('  Mean revisit time (covered pts): %.2f hr\n', mean(revisit_iss));
    fprintf('  Max revisit time (covered pts):  %.2f hr\n', max(revisit_iss));
end
fprintf('  Max passes at any point:   %d\n\n', max(cov_iss.n_passes(:)));

% Find coverage at Goddard's nearest grid point
[~, lat_idx] = min(abs(cov_iss.lat_vec - goddard_lat));
[~, lon_idx] = min(abs(cov_iss.lon_vec - goddard_lon));
iss_goddard_cov  = cov_iss.coverage_frac(lat_idx, lon_idx);
iss_goddard_rev  = cov_iss.revisit_mean_hr(lat_idx, lon_idx);
iss_goddard_pass = cov_iss.n_passes(lat_idx, lon_idx);
fprintf('  Goddard grid point (lat=%.0f, lon=%.0f):\n', ...
    cov_iss.lat_vec(lat_idx), cov_iss.lon_vec(lon_idx));
fprintf('    Coverage fraction: %.4f\n', iss_goddard_cov);
fprintf('    Passes over 1 day: %d\n', iss_goddard_pass);
if ~isnan(iss_goddard_rev)
    fprintf('    Mean revisit:      %.2f hr\n\n', iss_goddard_rev);
else
    fprintf('    Mean revisit:      N/A (fewer than 2 passes)\n\n');
end

fig_iss_cov = plotCoverage(cov_iss, 'Quantity', 'coverage', ...
    'GroundStations', goddard_gs);
saveas(fig_iss_cov, fullfile(OUT_DIR, 'iss_coverage_frac.png'));

fig_iss_rev = plotCoverage(cov_iss, 'Quantity', 'revisit', ...
    'GroundStations', goddard_gs);
saveas(fig_iss_rev, fullfile(OUT_DIR, 'iss_revisit_time.png'));

fprintf('  Saved: iss_coverage_frac.png, iss_revisit_time.png\n\n');

fprintf('=======================================================\n');
fprintf(' All outputs saved to: %s\n', OUT_DIR);
fprintf('=======================================================\n');
