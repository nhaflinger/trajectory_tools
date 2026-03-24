% example_dawn.m
% Reproduction of NASA's Dawn mission (2007-2018) using low-thrust spiral tools.
%
% Dawn trajectory:
%   Launch:          27 Sep 2007  (Delta II 7925H-9.5; provided Earth escape)
%   Mars flyby:      17 Feb 2009  (gravity assist ~0.38 km/s free ΔV)
%   Vesta arrival:   16 Jul 2011
%   Vesta departure: 05 Sep 2012
%   Ceres arrival:   06 Mar 2015  (end of prime mission)
%
% Propulsion: NSTAR ion engine — Isp ~3100 s, max thrust ~92 mN
% Spacecraft: wet mass ~1217.7 kg, xenon propellant ~425 kg
%
% Model:
%   The Delta II 7925H upper stage delivered Dawn to an Earth-escape
%   trajectory — the ion engines were NOT used for Earth escape.
%   Mars flyby is modelled as a patched-conic gravity assist.
%   All heliocentric cruise legs and the Vesta/Ceres spiral insertions
%   and departures are modelled using Edelbaum low-thrust spirals.

bodies = constants();
g0     = 9.80665e-3;   % km/s^2

outDir = fullfile(fileparts(mfilename('fullpath')), 'output_dawn');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('\n=================================================================\n');
fprintf('  NASA Dawn Mission Reproduction\n');
fprintf('  Launch: 27 Sep 2007  |  Ceres arrival: 06 Mar 2015\n');
fprintf('=================================================================\n\n');

%% ---- Mission timeline (Julian Dates) -----------------------------------
jd_launch        = julianDate(2007,  9, 27);
jd_mars_flyby    = julianDate(2009,  2, 17);
jd_vesta_arr     = julianDate(2011,  7, 16);
jd_vesta_dep     = julianDate(2012,  9,  5);
jd_ceres_arr     = julianDate(2015,  3,  6);

tof_launch_mars  = jd_mars_flyby - jd_launch;         % ~509 days
tof_mars_vesta   = jd_vesta_arr  - jd_mars_flyby;     % ~878 days
tof_vesta_ops    = jd_vesta_dep  - jd_vesta_arr;      % ~416 days (ops, not cruise)
tof_vesta_ceres  = jd_ceres_arr  - jd_vesta_dep;      % ~912 days

fprintf('Mission segments:\n');
fprintf('  Launch -> Mars flyby:    %.0f days (%.1f yr)\n', tof_launch_mars,  tof_launch_mars/365.25);
fprintf('  Mars flyby -> Vesta:     %.0f days (%.1f yr)\n', tof_mars_vesta,   tof_mars_vesta/365.25);
fprintf('  Vesta operations:        %.0f days (%.1f yr)\n', tof_vesta_ops,    tof_vesta_ops/365.25);
fprintf('  Vesta -> Ceres:          %.0f days (%.1f yr)\n', tof_vesta_ceres,  tof_vesta_ceres/365.25);
fprintf('  Total mission:           %.0f days (%.1f yr)\n\n', ...
    jd_ceres_arr - jd_launch, (jd_ceres_arr - jd_launch)/365.25);

%% ---- Spacecraft & engine parameters ------------------------------------
wet_mass   = 1217.7;   % kg  (Dawn at launch)
xenon_mass = 425.0;    % kg  (actual xenon loaded)
dry_mass   = wet_mass - xenon_mass;   % ~792.7 kg

isp_dawn   = 3100;     % s   (NSTAR Isp)
thrust_N   = 0.092;    % N   (NSTAR max thrust ~92 mN)

fprintf('Spacecraft parameters:\n');
fprintf('  Wet mass at launch:   %.1f kg\n', wet_mass);
fprintf('  Xenon propellant:     %.1f kg\n', xenon_mass);
fprintf('  Dry mass (estimate):  %.1f kg\n', dry_mass);
fprintf('  NSTAR Isp:            %d s\n', isp_dawn);
fprintf('  NSTAR thrust:         %.0f mN\n', thrust_N * 1000);
fprintf('\n');

%% ---- SOI / parking orbit radii -----------------------------------------
muSun       = bodies.Sun.mu;
r_SOI_earth = bodies.Earth.a * (bodies.Earth.mu / muSun)^(2/5);
r_SOI_vesta = bodies.Vesta.a * (bodies.Vesta.mu / muSun)^(2/5);
r_SOI_ceres = bodies.Ceres.a * (bodies.Ceres.mu / muSun)^(2/5);

% ---- Delta II escape ΔV -------------------------------------------------
% Dawn was delivered to Earth escape by the Delta II 7925H upper stage.
% Compute the Edelbaum ΔV from the launch parking orbit to Earth's SOI —
% this is what the launch vehicle provided, not the ion engine.
alt_earth_park = 185;   % km  (Dawn's parking orbit, Cape Canaveral)
r_park_earth   = bodies.Earth.radius + alt_earth_park;
lv_esc_opts.thrustN = 0;   % analytical only
lv_esc = lowThrustSpiral(bodies.Earth, r_park_earth, r_SOI_earth, lv_esc_opts);
lv_escape_dv = lv_esc.deltaV;   % km/s  (LEO circular -> Earth SOI Edelbaum)

fprintf('Launch vehicle Earth escape:\n');
fprintf('  Parking orbit alt: %d km\n', alt_earth_park);
fprintf('  Earth SOI:         %.0f km  (%.4f AU)\n', r_SOI_earth, ...
    r_SOI_earth / bodies.Constants.AU);
fprintf('  ΔV (LV, Edelbaum): %.3f km/s\n\n', lv_escape_dv);

% Dawn's operational orbits at Vesta & Ceres
% Vesta: LAMO ~210 km altitude (r ~ 473 km from centre)
% Ceres: LAMO ~375 km altitude (r ~ 845 km from centre)
alt_vesta_lamo = 210;    % km
alt_ceres_lamo = 375;    % km
r_vesta_lamo = bodies.Vesta.radius + alt_vesta_lamo;
r_ceres_lamo = bodies.Ceres.radius + alt_ceres_lamo;

fprintf('Body radii / SOI:\n');
fprintf('  Vesta radius: %.1f km  |  SOI: %.0f km  |  LAMO alt: %d km\n', ...
    bodies.Vesta.radius, r_SOI_vesta, alt_vesta_lamo);
fprintf('  Ceres radius: %.1f km  |  SOI: %.0f km  |  LAMO alt: %d km\n\n', ...
    bodies.Ceres.radius, r_SOI_ceres, alt_ceres_lamo);

% Base thruster options — update .wetMass per phase
dawnBase.thrustN        = 0.092;
dawnBase.isp            = 3100;
dawnBase.nStepsPerOrbit = 50;
dawnBase.nOutputPoints  = 500;

%% ---- Leg 1: Earth escape -> Mars flyby (heliocentric cruise) -----------
% Delta II provides Earth escape; ion engine operates during cruise to Mars.
% Model: heliocentric Edelbaum from Earth's actual position at launch
%        to Mars' actual position at flyby.
fprintf('--- Leg 1: Heliocentric cruise  Earth -> Mars flyby ---\n');

r1_earth = norm(orbitalState(bodies.Earth, jd_launch));
r2_mars  = norm(orbitalState(bodies.Mars,  jd_mars_flyby));

dawnBase.wetMass = wet_mass;
leg1 = lowThrustSpiral(bodies.Sun, r1_earth, r2_mars, dawnBase);

m_at_mars = wet_mass;
if ~isnan(leg1.finalMass), m_at_mars = leg1.finalMass; end

fprintf('  Earth heliocentric dist:  %.4f AU\n', r1_earth / bodies.Constants.AU);
fprintf('  Mars  heliocentric dist:  %.4f AU\n', r2_mars  / bodies.Constants.AU);
fprintf('  Edelbaum ΔV:    %.3f km/s\n', leg1.deltaV);
fprintf('  Propellant:     %.1f kg\n',   wet_mass - m_at_mars);
fprintf('  TOF (thruster): %.0f days\n\n', leg1.tofDays);

%% ---- Mars gravity assist -----------------------------------------------
% The flyby was used to raise Dawn's perihelion and change the orbital plane
% for the Vesta trajectory.  The "free" ΔV equivalent is roughly 0.38 km/s.
fprintf('--- Mars gravity assist (17 Feb 2009) ---\n');
fprintf('  Equivalent free ΔV:  ~0.38 km/s  (raises apo, changes inclination)\n');
fprintf('  Mass after flyby:    %.1f kg  (no propellant used)\n\n', m_at_mars);

%% ---- Leg 2: Mars flyby -> Vesta insertion ------------------------------
fprintf('--- Leg 2: Heliocentric cruise  Mars -> Vesta + insertion ---\n');

r3_vesta_helio = norm(orbitalState(bodies.Vesta, jd_vesta_arr));
r2_mars_dep    = norm(orbitalState(bodies.Mars,  jd_mars_flyby));

% Heliocentric: Mars heliocentric distance at flyby -> Vesta heliocentric distance
dawnBase.wetMass = m_at_mars;
leg2_helio = lowThrustSpiral(bodies.Sun, r2_mars_dep, r3_vesta_helio, dawnBase);
m_after_cruise2 = m_at_mars;
if ~isnan(leg2_helio.finalMass), m_after_cruise2 = leg2_helio.finalMass; end

% Vesta insertion: SOI -> LAMO
dawnBase.wetMass = m_after_cruise2;
leg2_insert = lowThrustSpiral(bodies.Vesta, r_SOI_vesta, r_vesta_lamo, dawnBase);
m_at_vesta = m_after_cruise2;
if ~isnan(leg2_insert.finalMass), m_at_vesta = leg2_insert.finalMass; end

fprintf('  Mars  heliocentric dist:  %.4f AU\n', r2_mars_dep    / bodies.Constants.AU);
fprintf('  Vesta heliocentric dist:  %.4f AU\n', r3_vesta_helio / bodies.Constants.AU);
fprintf('  Heliocentric ΔV:  %.3f km/s  (propellant: %.1f kg)\n', ...
    leg2_helio.deltaV, m_at_mars - m_after_cruise2);
fprintf('  Vesta insertion ΔV:  %.3f km/s  SOI -> LAMO %.0f km  (propellant: %.1f kg)\n', ...
    leg2_insert.deltaV, alt_vesta_lamo, m_after_cruise2 - m_at_vesta);
fprintf('  Mass at Vesta LAMO:  %.1f kg\n\n', m_at_vesta);

%% ---- Leg 3: Vesta departure -> Ceres insertion -------------------------
fprintf('--- Leg 3: Vesta departure + heliocentric cruise + Ceres insertion ---\n');

% Vesta departure: LAMO -> SOI
dawnBase.wetMass = m_at_vesta;
leg3_escape = lowThrustSpiral(bodies.Vesta, r_vesta_lamo, r_SOI_vesta, dawnBase);
m_after_vesta_esc = m_at_vesta;
if ~isnan(leg3_escape.finalMass), m_after_vesta_esc = leg3_escape.finalMass; end

% Heliocentric: Vesta -> Ceres
r4_vesta_dep  = norm(orbitalState(bodies.Vesta, jd_vesta_dep));
r5_ceres_arr  = norm(orbitalState(bodies.Ceres, jd_ceres_arr));

dawnBase.wetMass = m_after_vesta_esc;
leg3_helio = lowThrustSpiral(bodies.Sun, r4_vesta_dep, r5_ceres_arr, dawnBase);
m_after_cruise3 = m_after_vesta_esc;
if ~isnan(leg3_helio.finalMass), m_after_cruise3 = leg3_helio.finalMass; end

% Ceres insertion: SOI -> LAMO
dawnBase.wetMass = m_after_cruise3;
leg3_insert = lowThrustSpiral(bodies.Ceres, r_SOI_ceres, r_ceres_lamo, dawnBase);
m_at_ceres = m_after_cruise3;
if ~isnan(leg3_insert.finalMass), m_at_ceres = leg3_insert.finalMass; end

fprintf('  Vesta escape ΔV:    %.3f km/s  LAMO -> SOI  (propellant: %.1f kg)\n', ...
    leg3_escape.deltaV, m_at_vesta - m_after_vesta_esc);
fprintf('  Vesta heliocentric dist:  %.4f AU\n', r4_vesta_dep / bodies.Constants.AU);
fprintf('  Ceres heliocentric dist:  %.4f AU\n', r5_ceres_arr / bodies.Constants.AU);
fprintf('  Heliocentric ΔV:    %.3f km/s  (propellant: %.1f kg)\n', ...
    leg3_helio.deltaV, m_after_vesta_esc - m_after_cruise3);
fprintf('  Ceres insertion ΔV: %.3f km/s  SOI -> LAMO %.0f km  (propellant: %.1f kg)\n', ...
    leg3_insert.deltaV, alt_ceres_lamo, m_after_cruise3 - m_at_ceres);
fprintf('  Mass at Ceres LAMO:  %.1f kg\n\n', m_at_ceres);

%% ---- Mission totals & comparison with actual ---------------------------
prop_model  = wet_mass - m_at_ceres;
prop_actual = xenon_mass;   % 425 kg loaded (not all consumed)

dv_leg1       = leg1.deltaV;
dv_leg2       = leg2_helio.deltaV + leg2_insert.deltaV;
dv_leg3       = leg3_escape.deltaV + leg3_helio.deltaV + leg3_insert.deltaV;
dv_total      = dv_leg1 + dv_leg2 + dv_leg3;

fprintf('=================================================================\n');
fprintf('  MISSION SUMMARY\n');
fprintf('=================================================================\n');
fprintf('  %-42s %10s\n', 'Segment', 'ΔV (km/s)');
fprintf('  %s\n', repmat('-', 1, 55));
fprintf('  %-42s %10.3f\n', 'Leg 1: Earth -> Mars (heliocentric)',      dv_leg1);
fprintf('  %-42s %10s\n',   '  [Mars gravity assist]',                  '~free');
fprintf('  %-42s %10.3f\n', 'Leg 2: Mars -> Vesta (cruise+insertion)',   dv_leg2);
fprintf('  %-42s %10.3f\n', 'Leg 3: Vesta escape+cruise+Ceres insert',  dv_leg3);
fprintf('  %s\n', repmat('-', 1, 55));
fprintf('  %-42s %10.3f\n', 'TOTAL mission ΔV',                         dv_total);
fprintf('\n');
fprintf('  %-35s  %8.1f kg\n', 'Propellant used (model):',      prop_model);
fprintf('  %-35s  %8.1f kg\n', 'Xenon loaded (actual Dawn):',   prop_actual);
fprintf('  %-35s  %8.1f kg\n', 'Final mass at Ceres (model):',  m_at_ceres);
fprintf('  %-35s  %8.1f kg\n', 'Dry mass (estimated):',         dry_mass);
fprintf('\n');

prop_frac_model  = prop_model  / wet_mass;
prop_frac_actual = prop_actual / wet_mass;
fprintf('  Propellant fraction (model):   %.1f%%\n', 100*prop_frac_model);
fprintf('  Propellant fraction (actual):  %.1f%%  (xenon loaded / wet mass)\n', 100*prop_frac_actual);
fprintf('\n');

fprintf('  Notes:\n');
fprintf('   - Edelbaum model assumes slow-spiral (optimal continuous thrust).\n');
fprintf('   - Actual Dawn thrust level varied with solar distance (power-limited).\n');
fprintf('   - Actual trajectory included coast arcs to conserve xenon.\n');
fprintf('   - Mars gravity assist replaced ~0.38 km/s of propulsive ΔV.\n');
fprintf('   - Model ΔV is an analytical lower bound; actual may differ ~10-20%%.\n');
fprintf('\n');

%% ---- ΔV vs. direct transfer (no Mars flyby) ----------------------------
fprintf('--- Comparison: Direct Earth -> Vesta (no Mars flyby) ---\n');

r_earth_jd_launch = norm(orbitalState(bodies.Earth, jd_launch));
r_vesta_jd_direct = norm(orbitalState(bodies.Vesta, jd_vesta_arr));

dawnBase.wetMass = wet_mass;
leg_direct = lowThrustSpiral(bodies.Sun, r_earth_jd_launch, r_vesta_jd_direct, dawnBase);
m_direct_at_vesta_soi = wet_mass * exp(-leg_direct.deltaV / (isp_dawn * g0));
dawnBase.wetMass = m_direct_at_vesta_soi;
leg_direct_insert = lowThrustSpiral(bodies.Vesta, r_SOI_vesta, r_vesta_lamo, dawnBase);

dv_direct = leg_direct.deltaV + leg_direct_insert.deltaV;

fprintf('  Direct Earth -> Vesta ΔV:  %.3f km/s\n', dv_direct);
fprintf('  Dawn trajectory ΔV (E->V): %.3f km/s\n', dv_leg1 + dv_leg2);
fprintf('  Mars flyby saved:          %.3f km/s\n\n', dv_direct - (dv_leg1 + dv_leg2));

%% ---- Trajectory plots --------------------------------------------------
fprintf('Generating trajectory plots...\n');

% Build skipped-phase structs (no trajectory, zero ΔV/prop)
ph_skip_dep.deltaV        = 0;   ph_skip_dep.tof = 0;
ph_skip_dep.tofDays       = 0;   ph_skip_dep.propellantMass = 0;
ph_skip_dep.finalMass     = wet_mass;
ph_skip_dep.trajectory    = struct();

ph_skip_arr_mars.deltaV       = 0;   ph_skip_arr_mars.tof = 0;
ph_skip_arr_mars.tofDays      = 0;   ph_skip_arr_mars.propellantMass = 0;
ph_skip_arr_mars.finalMass    = m_at_mars;
ph_skip_arr_mars.trajectory   = struct();

% --- Plot 1: Earth -> Mars (Leg 1, heliocentric only) --------------------
res_leg1.deltaV         = leg1.deltaV;
res_leg1.propellantMass = wet_mass - m_at_mars;
res_leg1.finalMass      = m_at_mars;
res_leg1.tof            = NaN;
res_leg1.tofDays        = NaN;
res_leg1.details.dvDeparture    = 0;
res_leg1.details.dvHeliocentric = leg1.deltaV;
res_leg1.details.dvArrival      = 0;
res_leg1.details.r1_helio       = r1_earth;
res_leg1.details.r2_helio       = r2_mars;
res_leg1.phases.departure       = ph_skip_dep;
res_leg1.phases.heliocentric    = leg1;
res_leg1.phases.arrival         = ph_skip_arr_mars;
fig1 = plotLowThrustInterplanetary(res_leg1, bodies.Earth, bodies.Mars, ...
    struct('lvEscapeDV', lv_escape_dv, 'lvParkAlt', alt_earth_park));
saveas(fig1, fullfile(outDir, 'dawn_leg1_earth_mars.png'));

% --- Plot 2: Mars -> Vesta (Leg 2, cruise + insertion) -------------------
ph_skip_dep2.deltaV       = 0;   ph_skip_dep2.tof = 0;
ph_skip_dep2.tofDays      = 0;   ph_skip_dep2.propellantMass = 0;
ph_skip_dep2.finalMass    = m_at_mars;
ph_skip_dep2.trajectory   = struct();

res_leg2.deltaV         = leg2_helio.deltaV + leg2_insert.deltaV;
res_leg2.propellantMass = m_at_mars - m_at_vesta;
res_leg2.finalMass      = m_at_vesta;
res_leg2.tof            = NaN;
res_leg2.tofDays        = NaN;
res_leg2.details.dvDeparture    = 0;
res_leg2.details.dvHeliocentric = leg2_helio.deltaV;
res_leg2.details.dvArrival      = leg2_insert.deltaV;
res_leg2.details.r1_helio       = r2_mars_dep;
res_leg2.details.r2_helio       = r3_vesta_helio;
res_leg2.phases.departure       = ph_skip_dep2;
res_leg2.phases.heliocentric    = leg2_helio;
res_leg2.phases.arrival         = leg2_insert;
fig2 = plotLowThrustInterplanetary(res_leg2, bodies.Mars, bodies.Vesta);
saveas(fig2, fullfile(outDir, 'dawn_leg2_mars_vesta.png'));

% --- Plot 3: Vesta -> Ceres (Leg 3, escape + cruise + insertion) ---------
res_leg3.deltaV         = leg3_escape.deltaV + leg3_helio.deltaV + leg3_insert.deltaV;
res_leg3.propellantMass = m_at_vesta - m_at_ceres;
res_leg3.finalMass      = m_at_ceres;
res_leg3.tof            = NaN;
res_leg3.tofDays        = NaN;
res_leg3.details.dvDeparture    = leg3_escape.deltaV;
res_leg3.details.dvHeliocentric = leg3_helio.deltaV;
res_leg3.details.dvArrival      = leg3_insert.deltaV;
res_leg3.details.r1_helio       = r4_vesta_dep;
res_leg3.details.r2_helio       = r5_ceres_arr;
res_leg3.phases.departure       = leg3_escape;
res_leg3.phases.heliocentric    = leg3_helio;
res_leg3.phases.arrival         = leg3_insert;
fig3 = plotLowThrustInterplanetary(res_leg3, bodies.Vesta, bodies.Ceres);
saveas(fig3, fullfile(outDir, 'dawn_leg3_vesta_ceres.png'));
