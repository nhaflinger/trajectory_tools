% example_low_thrust_interplanetary.m
% Phase 2 low-thrust tools: interplanetary transfer budget + pork-chop.
%
% Computes a low-thrust Earth->Mars mission where the launch vehicle upper
% stage provides Earth escape (skipDepartureSpiral = true).  The ion engine
% performs only the heliocentric cruise and Mars capture spiral.
% Then sweeps the 2026-2028 launch window as a low-thrust pork-chop.

bodies = constants();

fprintf('\n=================================================================\n');
fprintf('  Low-Thrust Interplanetary: Earth -> Mars\n');
fprintf('=================================================================\n\n');

%% ---- Mission parameters --------------------------------------------
departJD = julianDate(2026, 1, 1);
tofDays  = 800;   % total mission TOF (days)

opts = struct( ...
    'thrustN',             0.100,  ...   % N  (100 mN ion engine)
    'isp',                 3000,   ...   % s
    'wetMass',             1000,   ...   % kg
    'departureAltitude',   200,    ...   % km (reference; not used — LV provides escape)
    'arrivalAltitude',     400,    ...   % km LMO
    'skipDepartureSpiral', true);        % LV upper stage delivers spacecraft to Earth escape

%% ---- Three-phase mission budget ------------------------------------
res = lowThrustInterplanetary(bodies.Earth, bodies.Mars, departJD, tofDays, opts);
plotLowThrustInterplanetary(res, bodies.Earth, bodies.Mars);

fprintf('Three-phase low-thrust Earth -> Mars\n');
fprintf('  Departure: %s   TOF budget: %d days\n', ...
    datestr(departJD - 1721058.5, 'dd mmm yyyy'), tofDays);
fprintf('\n');

fprintf('  Phase 1 — Earth escape (provided by launch vehicle upper stage)\n');
fprintf('    ΔV:    %.3f km/s  [skipped — LV provides escape]\n', res.details.dvDeparture);
fprintf('    Prop:  %.1f kg\n', res.details.propDeparture);
fprintf('    TOF:   %.0f days\n', res.details.tofDaysDeparture);

fprintf('  Phase 2 — Heliocentric cruise (Earth orbit -> Mars orbit)\n');
fprintf('    ΔV:    %.3f km/s\n', res.details.dvHeliocentric);
fprintf('    Prop:  %.1f kg\n', res.details.propHeliocentric);
fprintf('    TOF:   %.0f days\n', res.details.tofDaysHeliocentric);

fprintf('  Phase 3 — Mars capture spiral (SOI -> LMO)\n');
fprintf('    ΔV:    %.3f km/s\n', res.details.dvArrival);
fprintf('    Prop:  %.1f kg\n', res.details.propArrival);
fprintf('    TOF:   %.0f days\n', res.details.tofDaysArrival);

fprintf('\n  TOTAL\n');
fprintf('    ΔV:          %.3f km/s\n', res.deltaV);
fprintf('    Propellant:  %.1f kg  (%.1f%% of %.0f kg wet mass)\n', ...
    res.propellantMass, 100*res.propellantMass/opts.wetMass, opts.wetMass);
fprintf('    Final mass:  %.1f kg\n', res.finalMass);
fprintf('    Total TOF:   %.0f days\n\n', res.tofDays);

%% ---- Impulsive comparison (same departure, same TOF) ---------------
impOpts = struct( ...
    'departureJD',       departJD, ...
    'tofDays',           tofDays, ...
    'departureAltitude', 200, ...
    'arrivalAltitude',   400, ...
    'arrivalApogeeAltitude', 400);
imp = patchedConicTransfer(bodies.Earth, bodies.Mars, impOpts);

g0 = 9.80665e-3;
imp_propFrac_chem = 1 - exp(-imp.deltaVBurns / (450 * g0));
lt_propFrac       = res.propellantMass / opts.wetMass;

fprintf('  Comparison (same departure, same TOF = %d days):\n', tofDays);
fprintf('  %-35s %10s %15s\n', 'Metric', 'Impulsive', 'Low-Thrust');
fprintf('  %s\n', repmat('-', 1, 62));
fprintf('  %-35s %10.3f %15.3f km/s\n', 'Total ΔV',      imp.deltaVBurns, res.deltaV);
fprintf('  %-35s %10.1f %15.1f %%\n',   'Propellant fraction (%%)', ...
    100*imp_propFrac_chem, 100*lt_propFrac);
fprintf('\n');
fprintf('  Key insight: low-thrust requires MORE ΔV (no Oberth effect)\n');
fprintf('  but FAR LESS propellant mass because Isp = %d s vs ~450 s chemical.\n\n', opts.isp);

%% ---- Isp trade: same ΔV, different Isp ----------------------------
fprintf('  Propellant fraction vs. Isp for %.2f km/s mission ΔV:\n', res.deltaV);
fprintf('  %-20s  %8s  %12s\n', 'Propulsion type', 'Isp (s)', 'Prop fraction');
fprintf('  %s\n', repmat('-', 1, 44));
ispList   = [450,  800, 1800, 3000, 5000];
ispNames  = {'Chemical (biprop)', 'Advanced chem', 'Hall thruster', ...
             'Ion engine', 'High-Isp ion'};
for k = 1:numel(ispList)
    pf = 1 - exp(-res.deltaV / (ispList(k) * g0));
    fprintf('  %-20s  %8d  %12.1f %%\n', ispNames{k}, ispList(k), 100*pf);
end
fprintf('\n');

%% ---- Low-thrust pork-chop: 2026-2028 launch window ----------------
fprintf('Generating low-thrust pork-chop plot...\n');

jdStart   = julianDate(2026,  1,  1);
jdEnd     = julianDate(2028, 12, 31);
depDates  = linspace(jdStart, jdEnd, 60);
tofRange  = linspace(400, 1400, 50);

pcOpts = struct( ...
    'isp',               opts.isp,   ...
    'wetMass',           opts.wetMass, ...
    'departureAltitude', 200,        ...
    'arrivalAltitude',   400,        ...
    'showImpulsive',     true);

porkChopLowThrust(bodies.Earth, bodies.Mars, depDates, tofRange, pcOpts);
