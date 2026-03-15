% example_lunar_south_pole.m
% Lunar transfer designed for operations at the lunar south pole.
%
% The spacecraft departs from a 28.5-degree LEO (Kennedy Space Center
% latitude) and captures into a high-apogee polar orbit. The perilune
% is positioned above the south pole to support surface access, while
% the high apolune reduces orbital maintenance costs and provides
% extended communications geometry over the pole.
%
% Capture orbit geometry (Moon-centered):
%   - Inclination: ~90 deg (polar, relative to Moon's equatorial plane)
%   - Perilune:    100 km  — low-altitude pass over south pole for landing
%   - Apolune:     5000 km — high apogee for reduced station-keeping

bodies = constants();

options = struct();
options.departureAltitude     = 200;   % km, LEO parking orbit altitude
options.arrivalAltitude       = 100;   % km, perilune altitude above south pole
options.arrivalApogeeAltitude = 5000;  % km, apolune altitude
options.departureInclination  = 28.5;  % deg, LEO inclination (KSC)
options.arrivalInclination    = 90;    % deg, polar capture orbit
options.arrivalArgOfPeriapsis = 90;    % deg; perilune above north pole, apolune at south pole

result = patchedConicTransfer(bodies.Earth, bodies.Moon, options);

%% Derived capture orbit parameters
Moon      = bodies.Moon;
r_peri    = Moon.radius + options.arrivalAltitude;
r_apo     = Moon.radius + options.arrivalApogeeAltitude;
a_cap     = (r_peri + r_apo) / 2;
e_cap     = (r_apo - r_peri) / (r_apo + r_peri);
T_cap_hr  = 2*pi * sqrt(a_cap^3 / Moon.mu) / 3600;
v_peri    = sqrt(2*Moon.mu/r_peri - Moon.mu/a_cap);  % km/s at perilune
v_apo     = sqrt(2*Moon.mu/r_apo  - Moon.mu/a_cap);  % km/s at apolune

%% Print mission summary
fprintf('Lunar South Pole Mission — Transfer Summary\n');
fprintf('============================================\n');
fprintf('  Total delta-V:            %.3f km/s\n', result.deltaV);
fprintf('  Time of flight:           %.2f days\n', result.tof / bodies.Constants.day);
fprintf('  TLI delta-V:              %.3f km/s\n', result.details.dvTLI);
fprintf('  Capture delta-V:          %.3f km/s\n', result.details.dvCapture);
fprintf('  Departure inclination:    %.1f deg (LEO)\n', options.departureInclination);
fprintf('  Arrival inclination:      %.1f deg (polar)\n', options.arrivalInclination);
fprintf('\nCapture Orbit Parameters\n');
fprintf('  Perilune altitude:        %.0f km  (south pole approach)\n', options.arrivalAltitude);
fprintf('  Apolune altitude:         %.0f km\n',  options.arrivalApogeeAltitude);
fprintf('  Semi-major axis:          %.0f km\n',  a_cap);
fprintf('  Eccentricity:             %.4f\n',     e_cap);
fprintf('  Orbital period:           %.2f hours\n', T_cap_hr);
fprintf('  Speed at perilune:        %.3f km/s\n', v_peri);
fprintf('  Speed at apolune:         %.3f km/s\n', v_apo);

%% Plot
plotPatchedConic(result, bodies.Earth, bodies.Moon);
