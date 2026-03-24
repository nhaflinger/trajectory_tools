%EXAMPLE_ORBIT_LIFETIME  Demonstrate orbitLifetime, cwPropagate, cwRendezvous,
%   and frozenOrbit functions.
%
%   Sections:
%     1. ISS-like LEO decay
%     2. Ballistic coefficient trade
%     3. Molniya lifetime
%     4. CW propagation demo (closed orbit verification)
%     5. CW rendezvous
%     6. Frozen orbit survey table
%     7. Frozen vs non-frozen propagation comparison

clear; close all; clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
OUT_DIR    = fullfile(SCRIPT_DIR, 'output_orbit_lifetime');
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];

fprintf('============================================\n');
fprintf('  Orbit Lifetime and Relative Motion Demo  \n');
fprintf('============================================\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 1. ISS-like LEO decay  (circular 410 km, 51.6 deg, CdAm=0.02 m^2/kg)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('--- 1. ISS-like LEO Decay (410 km circular, CdAm=0.02) ---\n');

orb_iss = earthOrbit('circular', 410, 51.6);

[res_iss, fig_iss] = orbitLifetime(orb_iss, ...
    'CdAm',       0.02, ...
    'MaxYears',   30,   ...
    'StepOrbits', 10,   ...
    'Plot',       true);

fprintf('  ISS lifetime: %.1f days (%.2f years)\n', ...
    res_iss.lifetime_days, res_iss.lifetime_years);
fprintf('  Reentered: %d\n', res_iss.reentered);

saveas(fig_iss, fullfile(OUT_DIR, 'lifetime_ISS_410km.png'));
fprintf('  Saved: lifetime_ISS_410km.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 2. Ballistic coefficient trade  (same 410 km orbit)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('--- 2. Ballistic Coefficient Trade (410 km, 5 CdAm values) ---\n');

CdAm_vec    = linspace(0.005, 0.05, 5);
lifetime_yr = zeros(1, numel(CdAm_vec));

for k = 1:numel(CdAm_vec)
    res_k = orbitLifetime(orb_iss, ...
        'CdAm',       CdAm_vec(k), ...
        'MaxYears',   30,           ...
        'StepOrbits', 10);
    lifetime_yr(k) = res_k.lifetime_years;
    fprintf('  CdAm = %.4f m^2/kg -> lifetime = %.2f yr\n', CdAm_vec(k), lifetime_yr(k));
end

fig2 = figure('Color', bgCol, 'Position', [100 100 800 450]);
ax2  = axes('Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4, 'Parent', fig2);
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

plot(ax2, CdAm_vec, lifetime_yr, 'o-', 'Color', [0.30 0.75 0.93], ...
    'LineWidth', 2.0, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.30 0.75 0.93], 'MarkerEdgeColor', 'none');

xlabel(ax2, 'Cd \cdot A/m  (m^2/kg)', 'Color', txtCol);
ylabel(ax2, 'Lifetime (years)',         'Color', txtCol);
title(ax2,  'Orbital Lifetime vs Ballistic Coefficient (410 km circular)', 'Color', txtCol);

saveas(fig2, fullfile(OUT_DIR, 'lifetime_CdAm_trade.png'));
fprintf('  Saved: lifetime_CdAm_trade.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 3. Molniya lifetime  (high eccentricity, rapid circularization then decay)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('--- 3. Molniya Lifetime (e=0.74, CdAm=0.01) ---\n');

orb_mol = earthOrbit('molniya', 0);

[res_mol, fig_mol] = orbitLifetime(orb_mol, ...
    'CdAm',       0.01, ...
    'MaxYears',   30,   ...
    'StepOrbits', 10,   ...
    'Plot',       true);

fprintf('  Molniya initial: a=%.1f km, e=%.3f, peri=%.1f km, apo=%.1f km\n', ...
    orb_mol.a, orb_mol.e, orb_mol.alt_peri, orb_mol.alt_apo);
fprintf('  Molniya lifetime: %.1f days (%.2f years)\n', ...
    res_mol.lifetime_days, res_mol.lifetime_years);

saveas(fig_mol, fullfile(OUT_DIR, 'lifetime_Molniya.png'));
fprintf('  Saved: lifetime_Molniya.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 4. CW propagation demo  (closed orbit, rho=5 km)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('--- 4. CW Propagation Demo (closed orbit, rho=5 km) ---\n');

orb_chief = earthOrbit('circular', 410, 51.6);
n_chief   = 2*pi / orb_chief.period;   % rad/s
rho       = 5.0;                        % km (radial semi-axis)

% CW closed-orbit initial conditions (no secular along-track drift)
%   dR(t) =  rho*cos(n*t)
%   dS(t) = -2*rho*sin(n*t)
%   dW(t) =  (rho/2)*cos(n*t)
%   Closed-orbit condition: dVS_0 = -2*n*dR_0
dr_lvlh0 = [rho; 0; rho/2];           % km, [R; S; W]
dv_cw0   = [0; -2*n_chief*rho; 0];    % km/s, CW rotating-frame

% Propagate 5 orbital periods with cwPropagate
n_orbits = 5;
T_chief  = orb_chief.period;
t_cw     = linspace(0, n_orbits * T_chief, 1000);

[r_cw, ~] = cwPropagate(dr_lvlh0, dv_cw0, n_chief, t_cw);

dR_cw = r_cw(1,:);
dS_cw = r_cw(2,:);
dW_cw = r_cw(3,:);

% Analytical closed-form
dR_anal = rho * cos(n_chief * t_cw);
dS_anal = -2*rho * sin(n_chief * t_cw);
dW_anal = (rho/2) * cos(n_chief * t_cw);

max_err_R = max(abs(dR_cw - dR_anal));
max_err_S = max(abs(dS_cw - dS_anal));
max_err_W = max(abs(dW_cw - dW_anal));

fprintf('  Comparing cwPropagate vs closed-form CW ellipse:\n');
fprintf('  Expected: dR ± %.2f km, dS ± %.2f km, dW ± %.2f km\n', rho, 2*rho, rho/2);
fprintf('  Max error dR: %.2e km\n', max_err_R);
fprintf('  Max error dS: %.2e km\n', max_err_S);
fprintf('  Max error dW: %.2e km\n', max_err_W);
fprintf('  (Errors should be near machine epsilon ~1e-14 km)\n');

fig4 = figure('Color', bgCol, 'Position', [100 100 1100 500]);

ax4a = subplot(1,2,1, 'Parent', fig4);
set(ax4a, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
hold(ax4a, 'on'); grid(ax4a, 'on'); box(ax4a, 'on'); axis(ax4a, 'equal');

plot(ax4a, dS_cw,   dR_cw,   '-',  'Color', [0.30 0.75 0.93], 'LineWidth', 2.5, ...
    'DisplayName', 'cwPropagate STM');
plot(ax4a, dS_anal, dR_anal, '--', 'Color', [0.95 0.85 0.25], 'LineWidth', 1.5, ...
    'DisplayName', 'Closed-form');
plot(ax4a, dS_cw(1), dR_cw(1), 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.90 0.40 0.40], 'MarkerEdgeColor', 'none', 'DisplayName', 't=0');
plot(ax4a, 0, 0, '+', 'MarkerSize', 12, 'Color', [0.80 0.80 0.80], ...
    'LineWidth', 1.5, 'DisplayName', 'Chief');

xlabel(ax4a, '\DeltaS  Along-track (km)', 'Color', txtCol);
ylabel(ax4a, '\DeltaR  Radial (km)',       'Color', txtCol);
title(ax4a, sprintf('CW Closed Orbit (R-S Plane, %d orbits)', n_orbits), 'Color', txtCol);
leg4a = legend(ax4a, 'Location', 'northeast');
set(leg4a, 'TextColor', txtCol, 'Color', bgCol+0.05, 'EdgeColor', [0.4 0.4 0.4]);

ax4b = subplot(1,2,2, 'Parent', fig4);
set(ax4b, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
hold(ax4b, 'on'); grid(ax4b, 'on'); box(ax4b, 'on');

t_min4 = t_cw / 60;
plot(ax4b, t_min4, dR_cw, '-',  'Color', [0.30 0.75 0.93], 'LineWidth', 1.5, 'DisplayName', '\DeltaR radial');
plot(ax4b, t_min4, dS_cw, '-',  'Color', [0.95 0.60 0.20], 'LineWidth', 1.5, 'DisplayName', '\DeltaS along-track');
plot(ax4b, t_min4, dW_cw, '--', 'Color', [0.80 0.50 0.90], 'LineWidth', 1.2, 'DisplayName', '\DeltaW cross-track');

xlabel(ax4b, 'Time (min)',            'Color', txtCol);
ylabel(ax4b, 'Relative position (km)', 'Color', txtCol);
title(ax4b, 'CW Components vs Time', 'Color', txtCol);
leg4b = legend(ax4b, 'Location', 'northeast');
set(leg4b, 'TextColor', txtCol, 'Color', bgCol+0.05, 'EdgeColor', [0.4 0.4 0.4]);

saveas(fig4, fullfile(OUT_DIR, 'cw_closed_orbit.png'));
fprintf('  Saved: cw_closed_orbit.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 5. CW rendezvous  (deputy 10 km ahead along-track, 1 km above)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('--- 5. CW Rendezvous (tf = 0.5 orbital periods) ---\n');

% Chief at ISS orbit
n_rend  = n_chief;          % rad/s
T_rend  = T_chief;          % s

% Deputy initial state: 10 km ahead in S, 1 km above in R
dr0_rend = [1.0; 10.0; 0.0];    % km, [R; S; W]
dv0_rend = [0.0;  0.0; 0.0];    % km/s, starting at rest in LVLH frame

% Target: dock at chief (origin)
drf_rend = [0.0; 0.0; 0.0];
dvf_rend = [0.0; 0.0; 0.0];

% Transfer time: 0.75 orbital periods
% Note: 0.5*T is a resonance (tau=pi zeros the W diagonal of Phi_rv — singular).
% 0.75*T (tau=3pi/2) is well-conditioned.
tf_rend = 0.75 * T_rend;    % s

res_rend = cwRendezvous(dr0_rend, dv0_rend, n_rend, tf_rend, drf_rend, dvf_rend);

% Plot: R-S plane showing trajectory
fig5 = figure('Color', bgCol, 'Position', [100 100 800 700]);
ax5  = axes('Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4, 'Parent', fig5);
hold(ax5, 'on'); grid(ax5, 'on'); box(ax5, 'on');

% Plot rendezvous trajectory (R-S plane)
plot(ax5, res_rend.r_hist(2,:), res_rend.r_hist(1,:), '-', ...
    'Color', [0.30 0.75 0.93], 'LineWidth', 2.5, 'DisplayName', 'Rendezvous trajectory');

% Mark initial position
plot(ax5, dr0_rend(2), dr0_rend(1), 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.95 0.85 0.25], 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Deputy initial (t=0)');

% Mark docking point (chief location)
plot(ax5, 0, 0, 'p', 'MarkerSize', 14, ...
    'MarkerFaceColor', [0.30 0.90 0.50], 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Chief / docking target');

% Annotate dv1 and dv2
quiver(ax5, dr0_rend(2), dr0_rend(1), ...
    res_rend.dv1(2)*500, res_rend.dv1(1)*500, 0, ...
    'Color', [0.95 0.50 0.20], 'LineWidth', 2.0, 'MaxHeadSize', 0.6, ...
    'DisplayName', sprintf('\\Deltav_1 = %.1f m/s', res_rend.dv1_mag));

quiver(ax5, drf_rend(2), drf_rend(1), ...
    res_rend.dv2(2)*500, res_rend.dv2(1)*500, 0, ...
    'Color', [0.90 0.30 0.60], 'LineWidth', 2.0, 'MaxHeadSize', 0.6, ...
    'DisplayName', sprintf('\\Deltav_2 = %.1f m/s', res_rend.dv2_mag));

xlabel(ax5, '\DeltaS  Along-track (km)', 'Color', txtCol);
ylabel(ax5, '\DeltaR  Radial (km)',       'Color', txtCol);
title(ax5, sprintf('CW Rendezvous  (t_f = 0.5 T = %.0f s)\nTotal \\DeltaV = %.1f m/s', ...
    tf_rend, res_rend.total_dv_m_s), 'Color', txtCol);
leg5 = legend(ax5, 'Location', 'best');
set(leg5, 'TextColor', txtCol, 'Color', bgCol+0.05, 'EdgeColor', [0.4 0.4 0.4]);

saveas(fig5, fullfile(OUT_DIR, 'cw_rendezvous.png'));
fprintf('  Saved: cw_rendezvous.png\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 6. Frozen orbit survey table
%% ════════════════════════════════════════════════════════════════════════════
fprintf('--- 6. Frozen Orbit Survey ---\n');

alts = [500, 600, 700, 800];
incs = [63.4, 75, 90, 98];

fprintf('\n  %-8s  %-8s  %-14s  %-12s  %-12s  %-12s\n', ...
    'Alt(km)', 'Inc(deg)', 'e_frozen', 'alt_peri(km)', 'alt_apo(km)', 'omega(deg)');
fprintf('  %s\n', repmat('-', 1, 72));

for k_a = 1:numel(alts)
    for k_i = 1:numel(incs)
        frz_k = frozenOrbit(alts(k_a), incs(k_i));
        fprintf('  %-8.0f  %-8.2f  %-14.4e  %-12.4f  %-12.4f  %-12.1f\n', ...
            alts(k_a), incs(k_i), frz_k.e, ...
            frz_k.alt_peri, frz_k.alt_apo, frz_k.omega);
    end
end
fprintf('\n');

%% ════════════════════════════════════════════════════════════════════════════
%% 7. Frozen vs non-frozen propagation comparison (600 km, 98 deg, 30 days)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('--- 7. Frozen vs Non-frozen Eccentricity (600 km, 98 deg, 30 days) ---\n');

% Compute frozen orbit elements
frz_600 = frozenOrbit(600, 98);

% Build frozen orbit struct for propagation
orb_frozen = earthOrbit('coe', frz_600.a, frz_600.e, frz_600.i, ...
    0, frz_600.omega, 0);

% Non-frozen: same altitude, but circular (e=0)
orb_circ = earthOrbit('circular', 600, 98);

% Propagate 30 days with J2+J3 method
dur_30 = 30 * 86400;   % s

fprintf('  Propagating frozen orbit (J3 method, 30 days)...\n');
res_frz = propagateOrbit(orb_frozen, dur_30, 'Method', 'j3', ...
    'StepSize', 600, 'OutputCOE', true);

fprintf('  Propagating circular orbit (J3 method, 30 days)...\n');
res_circ = propagateOrbit(orb_circ, dur_30, 'Method', 'j3', ...
    'StepSize', 600, 'OutputCOE', true);

t_days_frz  = res_frz.t  / 86400;
t_days_circ = res_circ.t / 86400;

e_frz  = res_frz.e;
e_circ = res_circ.e;

fprintf('  Frozen orbit e range: [%.4e, %.4e]  (delta = %.4e)\n', ...
    min(e_frz), max(e_frz), max(e_frz) - min(e_frz));
fprintf('  Circular orbit e range: [%.4e, %.4e]  (delta = %.4e)\n', ...
    min(e_circ), max(e_circ), max(e_circ) - min(e_circ));

fig7 = figure('Color', bgCol, 'Position', [100 100 1000 500]);
ax7  = axes('Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
    'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4, 'Parent', fig7);
hold(ax7, 'on'); grid(ax7, 'on'); box(ax7, 'on');

plot(ax7, t_days_frz,  e_frz,  '-',  'Color', [0.30 0.75 0.93], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Frozen (e_0=%.4e, \\omega_0=90°)', frz_600.e));
plot(ax7, t_days_circ, e_circ, '--', 'Color', [0.95 0.60 0.20], 'LineWidth', 2.0, ...
    'DisplayName', 'Circular (e_0=0)');

xlabel(ax7, 'Time (days)', 'Color', txtCol);
ylabel(ax7, 'Eccentricity', 'Color', txtCol);
title(ax7, 'Frozen vs Non-frozen Orbit: J3 Eccentricity Evolution (600 km, 98 deg, 30 days)', ...
    'Color', txtCol);
leg7 = legend(ax7, 'Location', 'best');
set(leg7, 'TextColor', txtCol, 'Color', bgCol+0.05, 'EdgeColor', [0.4 0.4 0.4]);

saveas(fig7, fullfile(OUT_DIR, 'frozen_vs_circular_eccentricity.png'));
fprintf('  Saved: frozen_vs_circular_eccentricity.png\n\n');

fprintf('============================================\n');
fprintf('  All examples complete.\n');
fprintf('  Output saved to: %s\n', OUT_DIR);
fprintf('============================================\n');
