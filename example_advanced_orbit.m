%EXAMPLE_ADVANCED_ORBIT  Demonstrate new Earth orbit tools functionality.
%
%   Covers: sunPosition, betaAngle, SSO LTAN, GEO station-keeping,
%   topocentricAzEl, lvlhFrame, wgs84Geodetic, dogLegTrade,
%   and extended propagation methods (J3/J4/SRP).

clear; close all; clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
OUT_DIR    = fullfile(SCRIPT_DIR, 'output_advanced_orbit');
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];

fprintf('========================================\n');
fprintf('  Advanced Earth Orbit Tools Examples\n');
fprintf('========================================\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 1. Sun Position at J2000 and Today
%% ════════════════════════════════════════════════════════════════════════════
fprintf('--- 1. Sun Position ---\n');

jd_j2000 = 2451545.0;
[r_hat_j2000, r_km_j2000] = sunPosition(jd_j2000);
sun_RA_j2000  = atan2d(r_hat_j2000(2), r_hat_j2000(1));
sun_Dec_j2000 = asind(r_hat_j2000(3));
fprintf('  J2000 (JD %.1f):\n', jd_j2000);
fprintf('    RA  = %.4f deg\n', sun_RA_j2000);
fprintf('    Dec = %.4f deg\n', sun_Dec_j2000);
fprintf('    |r| = %.0f km (%.4f AU)\n', norm(r_km_j2000), norm(r_km_j2000)/149597870.7);

% Today: 2026-03-23
jd_today = julianDate(2026, 3, 23);
[r_hat_today, r_km_today] = sunPosition(jd_today);
sun_RA_today  = atan2d(r_hat_today(2), r_hat_today(1));
sun_Dec_today = asind(r_hat_today(3));
fprintf('\n  Today (2026-03-23, JD %.1f):\n', jd_today);
fprintf('    RA  = %.4f deg\n', sun_RA_today);
fprintf('    Dec = %.4f deg\n', sun_Dec_today);
fprintf('    |r| = %.0f km (%.4f AU)\n', norm(r_km_today), norm(r_km_today)/149597870.7);

%% ════════════════════════════════════════════════════════════════════════════
%% 2. Beta Angle for ISS
%% ════════════════════════════════════════════════════════════════════════════
fprintf('\n--- 2. Beta Angle: ISS (410 km, 51.6 deg) ---\n');

orb_iss = earthOrbit('circular', 410, 51.6);
beta_iss = betaAngle(orb_iss, 365);

fig2 = plotBetaAngle(beta_iss);
saveas(fig2, fullfile(OUT_DIR, 'beta_angle_ISS.png'));
fprintf('  Saved: beta_angle_ISS.png\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 3. Beta Angle for SSO with LTAN = 10:30
%% ════════════════════════════════════════════════════════════════════════════
fprintf('\n--- 3. Beta Angle: SSO 550 km, LTAN=10:30 ---\n');

epoch_jd_sso = julianDate(2024, 1, 1);
orb_sso = earthOrbit('sso', 550, 0, 90, 10.5, epoch_jd_sso);
fprintf('  SSO inclination: %.4f deg\n', orb_sso.i);
fprintf('  LTAN            : %.2f hr\n', orb_sso.ltan_hrs);
fprintf('  RAAN from LTAN  : %.2f deg\n', orb_sso.sun_sync_RAAN);

beta_sso = betaAngle(orb_sso, 365);
fprintf('  Beta range: %.2f to %.2f deg (shows slow, bounded variation)\n', ...
    beta_sso.beta_min_deg, beta_sso.beta_max_deg);

fig3 = plotBetaAngle(beta_sso);
saveas(fig3, fullfile(OUT_DIR, 'beta_angle_SSO_LTAN1030.png'));
fprintf('  Saved: beta_angle_SSO_LTAN1030.png\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 4. GEO Station-Keeping Budgets
%% ════════════════════════════════════════════════════════════════════════════
fprintf('\n--- 4. GEO Station-Keeping Budgets ---\n');

slots = [0, 90, -120];   % 0 deg E (Atlantic), 90 deg E (Asia), 240 deg E (120 W Americas)
slot_names = {'0 deg E (Atlantic)', '90 deg E (Asia)', '240 deg E (120 W, Americas)'};

fprintf('\n  %-28s  %10s  %10s  %10s\n', 'Slot', 'N-S m/s/yr', 'E-W m/s/yr', 'Total m/s/yr');
fprintf('  %s\n', repmat('-', 1, 65));

for k = 1:3
    bgt = geoStationKeeping(slots(k));
    fprintf('  %-28s  %10.2f  %10.2f  %10.2f\n', ...
        slot_names{k}, bgt.dv_ns_m_s, bgt.dv_ew_m_s, bgt.dv_total_m_s);
end

%% ════════════════════════════════════════════════════════════════════════════
%% 5. Topocentric Az/El from Kennedy Space Center
%% ════════════════════════════════════════════════════════════════════════════
fprintf('\n--- 5. Topocentric Az/El from KSC ---\n');

% KSC coordinates
ksc_lat =  28.524;   % deg N
ksc_lon = -80.651;   % deg E
ksc_alt =    0.0;    % km

% Propagate ISS for 3 orbits
orb_iss2 = earthOrbit('circular', 410, 51.6);
T_iss    = orb_iss2.period;   % seconds
res_iss  = propagateOrbit(orb_iss2, 3 * T_iss, 'Method', 'j2', 'StepSize', 30);

% Compute JD at each time step
jd_iss = orb_iss2.epoch_jd + res_iss.t / 86400;   % JD vector

[az, el, rng] = topocentricAzEl(ksc_lat, ksc_lon, ksc_alt, res_iss.r_eci, jd_iss);

% Find pass segments (el > 10 deg)
el_min = 10;   % deg
in_view = el > el_min;
t_hr    = res_iss.t / 3600;

fprintf('  Observer: KSC (lat=%.3f, lon=%.3f)\n', ksc_lat, ksc_lon);
fprintf('  Min elevation threshold: %.0f deg\n', el_min);

% Find pass start/end indices
transitions = diff([0, in_view, 0]);
pass_starts = find(transitions ==  1);
pass_ends   = find(transitions == -1) - 1;

fprintf('  Passes above %.0f deg: %d\n', el_min, numel(pass_starts));
for k = 1:numel(pass_starts)
    idx_range = pass_starts(k):pass_ends(k);
    [max_el_k, i_max] = max(el(idx_range));
    t_start = t_hr(pass_starts(k));
    t_end   = t_hr(pass_ends(k));
    az_at_max = az(idx_range(i_max));
    fprintf('  Pass %d: t=[%.2f,%.2f] hr, max el=%.1f deg at az=%.1f deg\n', ...
        k, t_start, t_end, max_el_k, az_at_max);
end

%% ════════════════════════════════════════════════════════════════════════════
%% 6. LVLH Relative Motion: Clohessy-Wiltshire Ellipse
%% ════════════════════════════════════════════════════════════════════════════
fprintf('\n--- 6. LVLH Relative Motion (CW ellipse, deputy vs ISS chief) ---\n');

% CW closed-orbit initial conditions produce a non-drifting 3D relative orbit:
%   dR(t) =  rho*cos(nt)           — radial oscillation
%   dS(t) = -2*rho*sin(nt)         — along-track (2x amplitude = 2:1 ellipse)
%   dW(t) = (rho/2)*cos(nt)        — cross-track (half amplitude, in phase with dR)
% Closed-orbit condition: dVS_0 = -2*n*dR_0 (no secular drift in S)

n_iss = 2*pi / orb_iss2.period;   % rad/s, ISS mean motion
rho   = 5.0;                       % km, formation semi-axis (radial)

% Initial relative state in rotating LVLH frame
dr_lvlh0 = [rho; 0; rho/2];            % [dR; dS; dW] km
dv_cw0   = [0; -2*n_iss*rho; 0];       % [dVR; dVS; dVW] km/s (CW velocities)

% Convert LVLH initial delta-state to ECI
%   v_deputy_ECI = v_chief_ECI + R_lvlh2eci * (dv_CW + omega x dr_LVLH)
%   where omega = n * W_hat = [0;0;n] in LVLH
DCM_eci2lvlh = lvlhFrame(orb_iss2.r_vec, orb_iss2.v_vec, 'dcm');
DCM_lvlh2eci = DCM_eci2lvlh';
omega_lvlh   = [0; 0; n_iss];
dr_eci0 = DCM_lvlh2eci * dr_lvlh0;
dv_eci0 = DCM_lvlh2eci * (dv_cw0 + cross(omega_lvlh, dr_lvlh0));

% Build deputy orbit from ECI state
orb_deputy = earthOrbit('eci', ...
    orb_iss2.r_vec + dr_eci0, ...
    orb_iss2.v_vec + dv_eci0, ...
    orb_iss2.epoch_jd);

% Propagate chief and deputy for 20 orbits using J2
% (J2 differential precession accumulates over many orbits, visibly separating
%  the numerical trajectory from the CW linear prediction)
n_orbits   = 20;
dur_norb   = n_orbits * orb_iss2.period;
res_chief  = propagateOrbit(orb_iss2,   dur_norb, 'Method', 'j2', 'StepSize', 60);
res_deputy = propagateOrbit(orb_deputy, dur_norb, 'Method', 'j2', 'StepSize', 60);

N_pts = min(size(res_chief.r_eci,2), size(res_deputy.r_eci,2));
dR = zeros(1,N_pts); dS = zeros(1,N_pts); dW = zeros(1,N_pts);

for k = 1:N_pts
    [dr_k, ~] = lvlhFrame(res_deputy.r_eci(:,k), res_deputy.v_eci(:,k), ...
                           res_chief.r_eci(:,k),  res_chief.v_eci(:,k), 'relative');
    dR(k) = dr_k(1); dS(k) = dr_k(2); dW(k) = dr_k(3);
end

t_s6   = res_chief.t(1:N_pts);
t_min6 = t_s6 / 60;

% CW analytical prediction (linear, no J2)
dR_cw =  rho * cos(n_iss * t_s6);
dS_cw = -2*rho * sin(n_iss * t_s6);
dW_cw = (rho/2) * cos(n_iss * t_s6);

fprintf('  Formation size rho = %.2f km, %d orbits (~%.1f hr)\n', rho, n_orbits, dur_norb/3600);
fprintf('  Expected (CW): dR ± %.2f km, dS ± %.2f km, dW ± %.2f km\n', rho, 2*rho, rho/2);
fprintf('  Numerical vs CW prediction (final-orbit max error):\n');
last_orb = t_s6 > (n_orbits-1)*orb_iss2.period;
fprintf('    dR: %.4f km\n', max(abs(dR(last_orb) - dR_cw(last_orb))));
fprintf('    dS: %.4f km\n', max(abs(dS(last_orb) - dS_cw(last_orb))));
fprintf('    dW: %.4f km\n', max(abs(dW(last_orb) - dW_cw(last_orb))));
fprintf('  (Deviation driven by nonlinear gravity + J2 differential precession)\n');

% Plot: left = R-S plane CW ellipse, right = time histories
fig6 = figure('Color', bgCol, 'Position', [100 100 1100 500]);

ax6a = subplot(1,2,1, 'Parent', fig6);
set(ax6a, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
hold(ax6a,'on'); grid(ax6a,'on'); box(ax6a,'on'); axis(ax6a,'equal');
plot(ax6a, dS,    dR,    '-',  'Color', [0.30 0.75 0.93], 'LineWidth', 2.5, 'DisplayName', 'J2 numerical');
plot(ax6a, dS_cw, dR_cw, '--', 'Color', [0.95 0.85 0.25], 'LineWidth', 1.5, 'DisplayName', 'CW analytical');
plot(ax6a, dS(1), dR(1), 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.90 0.40 0.40], 'MarkerEdgeColor', 'none', 'DisplayName', 't = 0');
plot(ax6a, 0, 0, '+', 'MarkerSize', 12, 'Color', [0.80 0.80 0.80], ...
    'LineWidth', 1.5, 'DisplayName', 'Chief (origin)');
xlabel(ax6a, '\DeltaS  Along-track (km)', 'Color', txtCol);
ylabel(ax6a, '\DeltaR  Radial (km)',       'Color', txtCol);
title(ax6a, sprintf('CW Relative Orbit  (R–S Plane, %d orbits)', n_orbits), 'Color', txtCol);
leg6a = legend(ax6a, 'Location', 'northeast');
set(leg6a, 'TextColor', txtCol, 'Color', bgCol+0.05, 'EdgeColor', [0.4 0.4 0.4]);

ax6b = subplot(1,2,2, 'Parent', fig6);
set(ax6b, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
hold(ax6b,'on'); grid(ax6b,'on'); box(ax6b,'on');
T_min = orb_iss2.period / 60;
for orb_k = 1:n_orbits
    xline(ax6b, orb_k*T_min, ':', 'Color', [0.35 0.35 0.40], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
plot(ax6b, t_min6, dR, '-',  'Color', [0.30 0.75 0.93], 'LineWidth', 1.5, 'DisplayName', '\DeltaR  radial (J2)');
plot(ax6b, t_min6, dS, '-',  'Color', [0.95 0.60 0.20], 'LineWidth', 1.5, 'DisplayName', '\DeltaS  along-track (J2)');
plot(ax6b, t_min6, dW, '--', 'Color', [0.80 0.50 0.90], 'LineWidth', 1.2, 'DisplayName', '\DeltaW  cross-track (J2)');
plot(ax6b, t_min6, dR_cw, ':', 'Color', [0.30 0.75 0.93]*0.6, 'LineWidth', 1.0, 'DisplayName', '\DeltaR  CW prediction');
plot(ax6b, t_min6, dS_cw, ':', 'Color', [0.95 0.60 0.20]*0.6, 'LineWidth', 1.0, 'DisplayName', '\DeltaS  CW prediction');
xlabel(ax6b, 'Time (min)', 'Color', txtCol);
ylabel(ax6b, 'Relative position (km)', 'Color', txtCol);
title(ax6b, sprintf('Deputy relative to Chief  (%d orbits)', n_orbits), 'Color', txtCol);
leg6b = legend(ax6b, 'Location', 'northeast');
set(leg6b, 'TextColor', txtCol, 'Color', bgCol+0.05, 'EdgeColor', [0.4 0.4 0.4]);

saveas(fig6, fullfile(OUT_DIR, 'lvlh_cw_ellipse.png'));
fprintf('  Saved: lvlh_cw_ellipse.png\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 7. WGS84 Geodetic vs Spherical Approximation
%% ════════════════════════════════════════════════════════════════════════════
fprintf('\n--- 7. WGS84 vs Spherical Geodetic ---\n');

% Test point: 45 deg N, 90 deg E, 500 km altitude
lat_test = 45.0;
lon_test = 90.0;
alt_test = 500.0;   % km

% WGS84 geodetic -> ECEF
r_ecef_wgs84 = wgs84Geodetic(lat_test, lon_test, alt_test, 'inverse');
fprintf('  Test point: lat=%.1f deg, lon=%.1f deg, alt=%.0f km\n', lat_test, lon_test, alt_test);
fprintf('  WGS84 ECEF: [%.4f, %.4f, %.4f] km\n', r_ecef_wgs84(1), r_ecef_wgs84(2), r_ecef_wgs84(3));

% Round-trip: ECEF -> WGS84 geodetic
[lat_rt, lon_rt, alt_rt] = wgs84Geodetic(r_ecef_wgs84);
fprintf('  Round-trip errors:\n');
fprintf('    lat error: %.2e deg\n', abs(lat_rt - lat_test));
fprintf('    lon error: %.2e deg\n', abs(lon_rt - lon_test));
fprintf('    alt error: %.2e km\n',  abs(alt_rt - alt_test));

% Compare with spherical approximation
R_E_sph = 6378.1363;   % km, spherical radius
r_mag_sph = R_E_sph + alt_test;
r_ecef_sph = r_mag_sph * [cosd(lat_test)*cosd(lon_test); cosd(lat_test)*sind(lon_test); sind(lat_test)];
dr_norm = norm(r_ecef_wgs84 - r_ecef_sph);
fprintf('  Spherical vs WGS84 ECEF position error: %.4f km\n', dr_norm);

% Show WGS84 vs spherical for a range of latitudes
lats = 0:10:90;
err_km = zeros(size(lats));
for k = 1:numel(lats)
    r_w = wgs84Geodetic(lats(k), 0, 0, 'inverse');
    r_s = R_E_sph * [cosd(lats(k)); 0; sind(lats(k))];
    err_km(k) = norm(r_w - r_s);
end
fprintf('  Max WGS84 vs spherical error (surface, 0-90 lat): %.4f km\n', max(err_km));

%% ════════════════════════════════════════════════════════════════════════════
%% 8. Dog-Leg Trade: KSC -> ISS Orbit
%% ════════════════════════════════════════════════════════════════════════════
fprintf('\n--- 8. Dog-Leg Trade: KSC to ISS orbit ---\n');

jd_launch = julianDate(2026, 3, 23, 12, 0, 0);
[dl_result, dl_fig] = dogLegTrade(ksc_lat, ksc_lon, 51.6, 260.0, jd_launch, ...
    'NDays', 1, 'InsertionAlt_km', 410, 'Plot', true);

if ~isempty(dl_fig)
    saveas(dl_fig, fullfile(OUT_DIR, 'dogleg_trade_KSC_ISS.png'));
    fprintf('  Saved: dogleg_trade_KSC_ISS.png\n');
end

%% ════════════════════════════════════════════════════════════════════════════
%% 9. Propagation Comparison: J2 vs J2+J3+J4 vs Drag
%% ════════════════════════════════════════════════════════════════════════════
fprintf('\n--- 9. Propagation Comparison: 30-day ISS altitude history ---\n');

dur_30days = 30 * 86400;   % seconds

fprintf('  Propagating j2 method...\n');
res_j2  = propagateOrbit(orb_iss2, dur_30days, 'Method', 'j2',   'StepSize', 600);

fprintf('  Propagating j4 method (J2)...\n');
res_j4  = propagateOrbit(orb_iss2, dur_30days, 'Method', 'j4',   'StepSize', 600);

fprintf('  Propagating drag method (J2+drag)...\n');
res_drg = propagateOrbit(orb_iss2, dur_30days, 'Method', 'drag', 'StepSize', 600, ...
                          'CdAm', 0.02);

t_days_30 = res_j2.t / 86400;

fig9 = figure('Color', bgCol, 'Position', [100 100 1000 500]);
ax9  = axes('Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
            'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
hold(ax9, 'on'); grid(ax9, 'on'); box(ax9, 'on');

plot(ax9, t_days_30, res_j2.alt,  '-',  'Color', [0.30 0.75 0.93], 'LineWidth', 1.5, 'DisplayName', 'J2 (analytical)');
plot(ax9, res_j4.t/86400, res_j4.alt, '-',  'Color', [0.95 0.60 0.20], 'LineWidth', 1.5, 'DisplayName', 'J2+J3+J4 (numerical)');
plot(ax9, res_drg.t/86400, res_drg.alt, '--', 'Color', [0.90 0.30 0.30], 'LineWidth', 1.5, 'DisplayName', 'J2+Drag (Cd*A/m=0.02 m^2/kg)');

xlabel(ax9, 'Time (days)', 'Color', txtCol);
ylabel(ax9, 'Altitude (km)', 'Color', txtCol);
title(ax9, 'ISS: 30-day altitude history by propagation method', 'Color', txtCol);
leg9 = legend(ax9, 'Location', 'southwest');
set(leg9, 'TextColor', txtCol, 'Color', bgCol + 0.05, 'EdgeColor', [0.4 0.4 0.4]);

% Print final altitudes
fprintf('  Final altitude (day 30):\n');
fprintf('    J2 analytical     : %.4f km\n', res_j2.alt(end));
fprintf('    J2 numerical: %.4f km\n', res_j4.alt(end));
fprintf('    J2+Drag           : %.4f km\n', res_drg.alt(end));
fprintf('  Total drag-induced decay over 30 days: %.4f km\n', res_drg.alt(1) - res_drg.alt(end));

saveas(fig9, fullfile(OUT_DIR, 'propagation_comparison_30day.png'));
fprintf('  Saved: propagation_comparison_30day.png\n');

fprintf('\n========================================\n');
fprintf('  All examples complete.\n');
fprintf('  Output saved to: %s\n', OUT_DIR);
fprintf('========================================\n');
