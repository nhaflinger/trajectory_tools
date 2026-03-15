% example_lunar_transfer.m
% Demonstrates a simple patched-conic lunar transfer using the provided
% constants and patchedConicTransfer function.

bodies = constants();

options = struct();
options.departureAltitude    = 200;   % km, Earth parking orbit altitude
options.arrivalAltitude      = 100;   % km, lunar orbit altitude
options.departureInclination = 28.5;  % deg, parking orbit inclination relative to lunar transfer plane
options.arrivalInclination   = 0;     % deg, lunar capture orbit inclination relative to arrival plane

result = patchedConicTransfer(bodies.Earth, bodies.Moon, options);

fprintf('Lunar Transfer (Earth->Moon)\n');
fprintf('  Total delta-V (approx): %.3f km/s\n', result.deltaV);
fprintf('  Time of flight (days): %.2f\n', result.tof / bodies.Constants.day);
fprintf('  TLI delta-V: %.3f km/s\n', result.details.dvTLI);
fprintf('  Capture delta-V: %.3f km/s\n', result.details.dvCapture);
fprintf('  Departure inclination: %.1f deg\n', result.details.departureInclination);
fprintf('  Arrival inclination:   %.1f deg\n', result.details.arrivalInclination);

% Plot results
plotPatchedConic(result, bodies.Earth, bodies.Moon);
