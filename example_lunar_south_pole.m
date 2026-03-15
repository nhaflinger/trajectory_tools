% example_lunar_south_pole.m
% Lunar transfer designed for operations at the lunar south pole.
%
% The spacecraft departs from a 28.5-degree LEO (Kennedy Space Center
% latitude) and captures into a high-apogee polar orbit. The perilune
% is positioned above the north pole (ω = 90°) so the descending pass
% approaches the south pole for landing, while the high apolune reduces
% orbital maintenance costs and provides extended communications geometry.
%
% Capture orbit geometry (Moon-centered):
%   - Inclination: ~90 deg (polar, relative to Moon's equatorial plane)
%   - Perilune:    100 km  — low-altitude pass; perilune above north pole
%   - Apolune:     5000 km — high apogee above south pole
%
% Two transfer modes are compared:
%   Direct     — single combined LOI + plane change at perilune
%   Bi-elliptic — three burns: equatorial capture, plane change at high
%                 apoapsis (Moon SOI), apoapsis trim to final orbit

bodies = constants();

%% Common options
options = struct();
options.departureAltitude     = 200;   % km, LEO parking orbit altitude
options.arrivalAltitude       = 100;   % km, perilune altitude
options.arrivalApogeeAltitude = 5000;  % km, apolune altitude
options.departureInclination  = 28.5;  % deg, LEO inclination (KSC)
options.arrivalInclination    = 90;    % deg, polar capture orbit
options.arrivalArgOfPeriapsis = 90;    % deg; perilune above north pole, apolune at south pole

%% Direct transfer (combined LOI + plane change)
options.transferMode = 'direct';
result_direct = patchedConicTransfer(bodies.Earth, bodies.Moon, options);

%% Bi-elliptic transfer (Moon SOI as intermediate apoapsis)
options.transferMode = 'biElliptic';
result_bi = patchedConicTransfer(bodies.Earth, bodies.Moon, options);

%% Derived capture orbit parameters
Moon      = bodies.Moon;
r_peri    = Moon.radius + options.arrivalAltitude;
r_apo     = Moon.radius + options.arrivalApogeeAltitude;
a_cap     = (r_peri + r_apo) / 2;
e_cap     = (r_apo - r_peri) / (r_apo + r_peri);
T_cap_hr  = 2*pi * sqrt(a_cap^3 / Moon.mu) / 3600;
v_peri    = sqrt(2*Moon.mu/r_peri - Moon.mu/a_cap);
v_apo     = sqrt(2*Moon.mu/r_apo  - Moon.mu/a_cap);

%% Print comparison
fprintf('\nLunar South Pole Mission — Transfer Comparison\n');
fprintf('================================================\n');
fprintf('%-30s  %10s  %10s\n', '', 'Direct', 'Bi-elliptic');
fprintf('%-30s  %10s  %10s\n', '', '------', '-----------');
fprintf('%-30s  %9.3f   %9.3f  km/s\n', 'TLI delta-V', ...
    result_direct.details.dvTLI, result_bi.details.dvTLI);
fprintf('%-30s  %9.3f   %9.3f  km/s\n', 'Total arrival delta-V', ...
    result_direct.details.dvCapture, result_bi.details.dvCapture);
fprintf('  %-28s  %9s   %9.3f  km/s\n', 'Burn 1 (equatorial LOI)', '--', result_bi.details.dvLOI);
fprintf('  %-28s  %9s   %9.3f  km/s\n', 'Burn 2 (plane change)', '--', result_bi.details.dvPlaneChange);
fprintf('  %-28s  %9s   %9.3f  km/s\n', 'Burn 3 (apoapsis trim)', '--', result_bi.details.dvApoapsisTrim);
fprintf('%-30s  %9.3f   %9.3f  km/s\n', 'TOTAL delta-V', ...
    result_direct.deltaV, result_bi.deltaV);
fprintf('%-30s  %9.2f   %9.2f  days\n', 'Time of flight', ...
    result_direct.tof / bodies.Constants.day, result_bi.tof / bodies.Constants.day);
fprintf('%-30s  %9.1f   %9.1f  km\n', 'Intermediate apoapsis alt', ...
    0, result_bi.details.rBiEllipticApoapsis - Moon.radius);
fprintf('%-30s  %+9.3f\n', 'Bi-elliptic savings', ...
    result_direct.deltaV - result_bi.deltaV);

fprintf('\nCapture Orbit Parameters\n');
fprintf('  Perilune altitude:        %.0f km\n',   options.arrivalAltitude);
fprintf('  Apolune altitude:         %.0f km\n',   options.arrivalApogeeAltitude);
fprintf('  Semi-major axis:          %.0f km\n',   a_cap);
fprintf('  Eccentricity:             %.4f\n',      e_cap);
fprintf('  Orbital period:           %.2f hours\n', T_cap_hr);
fprintf('  Speed at perilune:        %.3f km/s\n', v_peri);
fprintf('  Speed at apolune:         %.3f km/s\n', v_apo);

%% Trajectory plots (direct result — overview and body-centric views)
plotPatchedConic(result_direct, bodies.Earth, bodies.Moon);

%% Delta-V comparison figure
figure('Name', 'Delta-V Comparison: Direct vs Bi-elliptic', 'NumberTitle', 'off', ...
       'Position', [100 100 620 420]);

labels  = {'TLI', 'LOI (equatorial)', 'Plane change', 'Apoapsis trim'};
dv_dir  = [result_direct.details.dvTLI,  result_direct.details.dvCapture, 0, 0];
dv_bi   = [result_bi.details.dvTLI, result_bi.details.dvLOI, ...
           result_bi.details.dvPlaneChange, result_bi.details.dvApoapsisTrim];

x = 1:4;
b1 = bar(x - 0.2, dv_dir, 0.35, 'FaceColor', [0.2 0.5 1.0]);
hold on;
b2 = bar(x + 0.2, dv_bi,  0.35, 'FaceColor', [0.9 0.45 0.1]);
hold off;

set(gca, 'XTick', x, 'XTickLabel', labels);
ylabel('\DeltaV (km/s)');
title('Burn Budget: Direct vs Bi-elliptic (90° polar capture)');
legend([b1 b2], ...
    sprintf('Direct  (total %.3f km/s)',      result_direct.deltaV), ...
    sprintf('Bi-elliptic  (total %.3f km/s)', result_bi.deltaV), ...
    'Location', 'northwest');
grid on;
text(0.98, 0.97, sprintf('Savings: %.3f km/s', result_direct.deltaV - result_bi.deltaV), ...
     'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.1 0.6 0.1]);
