%% example_maneuvers.m
% Demonstrates hohmannTransfer, planeChange, phasingManeuver, deorbitBurn,
% and launchWindow / plotLaunchWindow for representative Earth orbit maneuvers.

SCRIPT_DIR = fileparts(mfilename('fullpath'));
OUT_DIR    = fullfile(SCRIPT_DIR, 'output_maneuvers');
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

%% ══════════════════════════════════════════════════════════════════════════════
%% 1. Hohmann Transfer: LEO to GEO
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  DEMO 1: Hohmann Transfer — LEO (300 km) to GEO (35786 km)\n');
fprintf('%s\n', repmat('=', 1, 70));

res_hoh = hohmannTransfer(300, 35786);

fprintf('\n  Summary:\n');
fprintf('    DV1      = %+.4f km/s  (%+.2f m/s)\n', res_hoh.dv1_km_s, res_hoh.dv1_km_s*1000);
fprintf('    DV2      = %+.4f km/s  (%+.2f m/s)\n', res_hoh.dv2_km_s, res_hoh.dv2_km_s*1000);
fprintf('    DV total = %.4f km/s  (%.2f m/s)\n', res_hoh.dv_total_km_s, res_hoh.dv_total_km_s*1000);
fprintf('    TOF      = %.2f h  (%.1f min)\n', res_hoh.tof_s/3600, res_hoh.tof_min);

fprintf('\n--- Bi-elliptic via 70000 km intermediate ---\n');
R_E = 6378.1363;
r_int = R_E + 70000;   % geocentric radius of intermediate apoapsis
res_be = hohmannTransfer(300, 35786, 'Bielliptic', r_int);

fprintf('\n  Summary:\n');
fprintf('    DV1      = %+.4f km/s\n', res_be.dv1_km_s);
fprintf('    DV2      = %+.4f km/s\n', res_be.dv2_km_s);
fprintf('    DV3      = %+.4f km/s\n', res_be.dv3_km_s);
fprintf('    DV total = %.4f km/s  (Hohmann was %.4f km/s)\n', ...
        res_be.dv_total_km_s, res_hoh.dv_total_km_s);
fprintf('    TOF      = %.2f h\n', res_be.tof_s/3600);

%% ══════════════════════════════════════════════════════════════════════════════
%% 2. ISS Reboost: 380 km → 410 km
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  DEMO 2: ISS Reboost — 380 km to 410 km\n');
fprintf('%s\n', repmat('=', 1, 70));

res_reboost = hohmannTransfer(380, 410);

fprintf('\n  Summary:\n');
fprintf('    DV1       = %+.4f km/s  (%+.3f m/s)\n', res_reboost.dv1_km_s, res_reboost.dv1_km_s*1000);
fprintf('    DV2       = %+.4f km/s  (%+.3f m/s)\n', res_reboost.dv2_km_s, res_reboost.dv2_km_s*1000);
fprintf('    DV total  = %.4f km/s  (%.3f m/s)  <-- small reboost burn\n', ...
        res_reboost.dv_total_km_s, res_reboost.dv_total_km_s*1000);
fprintf('    TOF       = %.2f min\n', res_reboost.tof_min);

%% ══════════════════════════════════════════════════════════════════════════════
%% 3. Plane Change: ISS-like orbit, 5 deg inclination change
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  DEMO 3: Plane Change — ISS-like orbit (410 km, 51.6 deg), Delta-i = 5 deg\n');
fprintf('%s\n', repmat('=', 1, 70));

iss = earthOrbit('circular', 410, 51.6);

fprintf('\n-- Pure plane change at 410 km --\n');
res_pc_pure = planeChange(iss, 5.0);

fprintf('\n  DV (pure plane change) = %.4f km/s  (%.2f m/s)\n', ...
        res_pc_pure.dv_km_s, res_pc_pure.dv_km_s*1000);

fprintf('\n-- Combined plane change + altitude raise to 420 km --\n');
res_pc_comb = planeChange(iss, 5.0, 'Type', 'combined', 'AltFinal_km', 420);

fprintf('\n  DV combined  = %.4f km/s  (%.2f m/s)\n', ...
        res_pc_comb.dv_km_s, res_pc_comb.dv_km_s*1000);
fprintf('  DV separate  = %.4f km/s  (%.2f m/s)\n', ...
        res_pc_comb.dv_separate_km_s, res_pc_comb.dv_separate_km_s*1000);
fprintf('  Savings      = %.4f km/s  (%.2f m/s)\n', ...
        res_pc_comb.savings_km_s, res_pc_comb.savings_km_s*1000);

%% ══════════════════════════════════════════════════════════════════════════════
%% 4. Phasing Maneuver: SSO constellation slot (45 deg catch-up, 5 phasing orbits)
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  DEMO 4: Phasing Maneuver — SSO 550 km, catch up 45 deg in 5 orbits\n');
fprintf('%s\n', repmat('=', 1, 70));

sso = earthOrbit('sso', 550);

res_phase = phasingManeuver(sso, 45.0, 5);

fprintf('\n  Summary:\n');
fprintf('    Phasing orbit alt range  : %.1f km to %.1f km\n', ...
        res_phase.alt_peri_phase_km, res_phase.alt_apo_phase_km);
fprintf('    Phasing period           : %.2f min\n', res_phase.period_phase_min);
fprintf('    DV per burn              : %.4f km/s  (%.3f m/s)\n', ...
        res_phase.dv_each_km_s, res_phase.dv_each_km_s*1000);
fprintf('    DV total                 : %.4f km/s  (%.3f m/s)\n', ...
        res_phase.dv_total_km_s, res_phase.dv_total_km_s*1000);
fprintf('    Time to complete         : %.2f hr  (%.1f min)\n', ...
        res_phase.time_to_complete_hr, res_phase.time_to_complete_s/60);

%% ══════════════════════════════════════════════════════════════════════════════
%% 5. Deorbit Burn: from ISS altitude (410 km)
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  DEMO 5: Deorbit Burn — from ISS altitude (410 km) to 80 km perigee\n');
fprintf('%s\n', repmat('=', 1, 70));

res_deorb = deorbitBurn(iss);  % default target perigee = 80 km

fprintf('\n  Summary:\n');
fprintf('    DV (retrograde)     = %.4f km/s  (%.2f m/s)\n', ...
        res_deorb.dv_km_s, res_deorb.dv_m_s);
fprintf('    TOF to 80 km peri   = %.2f min  (%.1f s)\n', ...
        res_deorb.tof_to_peri_min, res_deorb.tof_to_peri_s);
fprintf('    25-year compliant   = %s\n', mat2str(res_deorb.compliant_25yr));

%% ══════════════════════════════════════════════════════════════════════════════
%% 6. Launch Windows: KSC to ISS orbit
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  DEMO 6: Launch Windows — KSC to ISS orbit (i=51.6 deg, RAAN=260 deg)\n');
fprintf('%s\n', repmat('=', 1, 70));

ksc_lat =  28.573;
ksc_lon = -80.649;
iss_inc  =  51.6;
iss_raan = 260.0;
j2000    =  2451545.0;

wins_ksc = launchWindow(ksc_lat, ksc_lon, iss_inc, iss_raan, j2000, ...
                        'NDays', 1, 'SearchStep_s', 30);

fig6 = plotLaunchWindow(wins_ksc, ksc_lat, ksc_lon, iss_inc, iss_raan, j2000, 1);
saveas(fig6, fullfile(OUT_DIR, 'launch_windows_ksc_iss.png'));
fprintf('\n  Plot saved: launch_windows_ksc_iss.png\n');

%% ══════════════════════════════════════════════════════════════════════════════
%% 7. Launch Windows: Vandenberg to SSO (550 km)
%% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  DEMO 7: Launch Windows — Vandenberg to SSO 550 km (i~97.4 deg, RAAN=0 deg)\n');
fprintf('%s\n', repmat('=', 1, 70));

vafb_lat =  34.632;
vafb_lon = -120.611;
sso_inc  =  sso.i;    % computed SSO inclination for 550 km
sso_raan =  0.0;

wins_vafb = launchWindow(vafb_lat, vafb_lon, sso_inc, sso_raan, j2000, ...
                         'NDays', 1, 'SearchStep_s', 30);

fig7 = plotLaunchWindow(wins_vafb, vafb_lat, vafb_lon, sso_inc, sso_raan, j2000, 1);
saveas(fig7, fullfile(OUT_DIR, 'launch_windows_vafb_sso.png'));
fprintf('\n  Plot saved: launch_windows_vafb_sso.png\n');

fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('  All outputs saved to: %s\n', OUT_DIR);
fprintf('%s\n\n', repmat('=', 1, 70));
