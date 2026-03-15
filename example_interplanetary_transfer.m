% example_interplanetary_transfer.m
% Demonstrates a simple patched-conic transfer between two planets.

bodies = constants();

options = struct();
options.departureAltitude    = 200;   % km, Earth parking orbit altitude
options.arrivalAltitude      = 400;   % km, Mars capture orbit altitude
options.departureInclination = 28.5;  % deg, parking orbit inclination relative to ecliptic transfer plane
options.arrivalInclination   = 0;     % deg, capture orbit inclination relative to arrival plane

% Example: Earth -> Mars
result = patchedConicTransfer(bodies.Earth, bodies.Mars, options);

fprintf('Interplanetary Transfer (Earth->Mars)\n');
fprintf('  Total delta-V (approx): %.3f km/s\n', result.deltaV);
fprintf('  Time of flight (days): %.2f\n', result.tof / bodies.Constants.day);
fprintf('  Required phase angle (deg): %.1f\n', result.phaseAngle);
fprintf('  Departure DV: %.3f km/s\n', result.details.dvDeparture);
fprintf('  Arrival DV: %.3f km/s\n', result.details.dvArrival);
fprintf('  Departure inclination: %.1f deg\n', result.details.departureInclination);
fprintf('  Arrival inclination:   %.1f deg\n', result.details.arrivalInclination);

% Pork chop plot (approximate, circular orbits)
% departure window: 2026-01-01 through 2032-12-31
jdStart = julianDate(2026, 1, 1);
jdEnd   = julianDate(2032, 12, 31);
tofDays = linspace(120, 350, 80)';

best = findBestLaunchDate(bodies.Earth, bodies.Mars, jdStart, jdEnd, tofDays);
fprintf('Best 2026-2032 launch: JD %.1f (~%s) for TOF %.1f days (DV %.2f km/s)\n', ...
    best.departureJD, datestr(best.departureJD - 1721058.5, 'yyyy-mm-dd'), best.tofDays, best.deltaV);

porkChopPlot(bodies.Earth, bodies.Mars, best.departJD, tofDays);

% Update transfer phase to match best departure date
[rE, ~] = orbitalStateCircular(bodies.Earth, best.departureJD);
[rM, ~] = orbitalStateCircular(bodies.Mars, best.departureJD);
phaseAngle = atan2d(rM(2)-rE(2), rM(1)-rE(1));
result.phaseAngle = phaseAngle;

% Store selected launch date so plot can label it
result.departureJD = best.departureJD;

% Re-plot trajectory with adjusted phase
plotPatchedConic(result, bodies.Earth, bodies.Mars);
