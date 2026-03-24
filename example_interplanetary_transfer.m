% example_interplanetary_transfer.m
% Earth -> Mars transfer using patched-conic with Lambert solver.
% Demonstrates: eccentric planet orbits, C3, TCM budget, Type I vs II.

bodies = constants();

outDir = fullfile(fileparts(mfilename('fullpath')), 'output_interplanetary_transfer');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

options = struct();
options.departureAltitude    = 200;   % km, Earth LEO parking orbit
options.arrivalAltitude      = 400;   % km, Mars capture orbit altitude
options.arrivalApogeeAltitude= 400;   % km, circular capture
options.departureInclination = 28.5;  % deg, KSC parking orbit vs ecliptic
options.arrivalInclination   = 0;     % deg, capture orbit vs arrival plane

% ---- Pork chop: find best launch window (2026-2032) ----
jdStart = julianDate(2026, 1, 1);
jdEnd   = julianDate(2032, 12, 31);
tofDays = linspace(120, 350, 80)';

best = findBestLaunchDate(bodies.Earth, bodies.Mars, jdStart, jdEnd, tofDays);
fprintf('Best 2026-2032 launch: %s  TOF = %.0f days  DV_helio = %.3f km/s\n', ...
    datestr(best.departureJD - 1721058.5, 'yyyy-mm-dd'), best.tofDays, best.deltaV);

porkChopPlot(bodies.Earth, bodies.Mars, best.departJD, tofDays);
saveas(gcf, fullfile(outDir, 'pork_chop_earth_mars.png'));

% ---- High-fidelity Lambert transfer at best date ----
options.departureJD  = best.departureJD;
options.tofDays      = best.tofDays;
options.transferType = 'type1';   % short-way (<180 deg)
result = patchedConicTransfer(bodies.Earth, bodies.Mars, options);
result.departureJD = best.departureJD;

arrJD  = best.departureJD + best.tofDays;
fprintf('\nEarth -> Mars  (%s)\n', datestr(best.departureJD - 1721058.5, 'yyyy-mm-dd'));
fprintf('  Transfer type    : %s\n',   options.transferType);
fprintf('  TOF              : %.1f days\n', result.tof/86400);
fprintf('  Departure DV     : %.3f km/s\n', result.details.dvDeparture);
fprintf('  Arrival DV       : %.3f km/s\n', result.details.dvArrival);
fprintf('  TCM budget       : %.0f m/s\n',  result.details.dvTCM*1000);
fprintf('  Total DV (burns) : %.3f km/s\n', result.deltaVBurns);
fprintf('  Total DV (+TCM)  : %.3f km/s\n', result.deltaV);
fprintf('  C3               : %.2f km^2/s^2\n', result.details.C3);
fprintf('  v_inf departure  : %.3f km/s\n', result.details.vInfDepart);
fprintf('  v_inf arrival    : %.3f km/s\n', result.details.vInfArrive);
fprintf('  Arrival date     : %s\n', datestr(arrJD - 1721058.5, 'yyyy-mm-dd'));

plotPatchedConic(result, bodies.Earth, bodies.Mars);
saveas(gcf, fullfile(outDir, 'patched_conic_earth_mars.png'));

% ---- Type II comparison ----
options.transferType = 'type2';
result2 = patchedConicTransfer(bodies.Earth, bodies.Mars, options);
fprintf('\nType II comparison (long-way):\n');
fprintf('  Total DV (+TCM)  : %.3f km/s\n', result2.deltaV);
fprintf('  C3               : %.2f km^2/s^2\n', result2.details.C3);
