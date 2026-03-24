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
Earth  = bodies.Earth;
Moon   = bodies.Moon;

outDir = fullfile(fileparts(mfilename('fullpath')), 'output_lunar_south_pole');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% Common options
options = struct();
options.departureAltitude     = 200;   % km, LEO parking orbit altitude
options.arrivalAltitude       = 100;   % km, perilune altitude
options.arrivalApogeeAltitude = 5000;  % km, apolune altitude
options.departureInclination  = 28.5;  % deg, LEO inclination (KSC)
options.arrivalInclination    = 90;    % deg, polar capture orbit
options.arrivalArgOfPeriapsis = 90;    % deg; perilune above north pole, apolune at south pole

%% Compute transfers
options.transferMode = 'direct';
result_direct = patchedConicTransfer(Earth, Moon, options);

options.transferMode = 'biElliptic';
result_bi = patchedConicTransfer(Earth, Moon, options);

%% Derived transfer orbit quantities (TLI is identical for both modes)
r0           = result_direct.details.r0;
a_xfr        = result_direct.details.transferSemiMajor;
e_xfr        = 1 - r0 / a_xfr;
T_xfr_days   = 2*pi * sqrt(a_xfr^3 / Earth.mu) / 86400;
v_park       = sqrt(Earth.mu / r0);
v_perigee    = sqrt(2*Earth.mu/r0 - Earth.mu/a_xfr);

%% Derived arrival quantities
r_cap          = result_direct.details.rParkArrive;
v_inf_dir      = result_direct.details.vInf;
v_inf_bi       = result_bi.details.vInf;
v_perilune_dir = sqrt(v_inf_dir^2 + 2*Moon.mu/r_cap);
v_perilune_bi  = sqrt(v_inf_bi^2  + 2*Moon.mu/r_cap);

%% TCM budget (2% of TLI ΔV, minimum 10 m/s)
dv_tcm_dir   = max(0.010, 0.02 * result_direct.details.dvTLI);
dv_tcm_bi    = max(0.010, 0.02 * result_bi.details.dvTLI);
dv_total_dir = result_direct.deltaV + dv_tcm_dir;
dv_total_bi  = result_bi.deltaV     + dv_tcm_bi;

%% Capture orbit parameters
r_peri   = Moon.radius + options.arrivalAltitude;
r_apo    = Moon.radius + options.arrivalApogeeAltitude;
a_cap    = (r_peri + r_apo) / 2;
e_cap    = (r_apo - r_peri) / (r_apo + r_peri);
T_cap_hr = 2*pi * sqrt(a_cap^3 / Moon.mu) / 3600;
v_peri   = sqrt(2*Moon.mu/r_peri - Moon.mu/a_cap);
v_apo    = sqrt(2*Moon.mu/r_apo  - Moon.mu/a_cap);

%% ---- Console output ----
fprintf('\nLunar South Pole Mission  (Earth -> Moon)\n');
fprintf('==========================================\n\n');

fprintf('Departure — Earth\n');
fprintf('  Parking orbit altitude:      %5.0f km\n',    options.departureAltitude);
fprintf('  Parking orbit inclination:   %5.1f deg\n',   options.departureInclination);
fprintf('  Parking orbit speed:         %7.3f km/s\n',  v_park);
fprintf('  TLI injection speed:         %7.3f km/s\n',  v_perigee);
fprintf('  TLI delta-V:                 %7.3f km/s\n',  result_direct.details.dvTLI);

fprintf('\nTransfer Orbit (identical for both modes)\n');
fprintf('  Semi-major axis:             %7.0f km\n',    a_xfr);
fprintf('  Eccentricity:                %7.4f\n',        e_xfr);
fprintf('  Full period:                 %7.2f days\n',  T_xfr_days);
fprintf('  Time of flight (half):       %7.2f days\n',  result_direct.tof / 86400);

fprintf('\nArrival — Moon\n');
fprintf('  %-34s  %9s  %9s\n', '',               'Direct', 'Bi-elliptic');
fprintf('  %-34s  %9s  %9s\n', '',               '------', '-----------');
fprintf('  %-34s  %7.3f    %7.3f  km/s\n', 'v∞ at Moon SOI', ...
    v_inf_dir, v_inf_bi);
fprintf('  %-34s  %7.3f    %7.3f  km/s\n', 'Approach speed at perilune', ...
    v_perilune_dir, v_perilune_bi);
fprintf('  %-34s  %5.0f\n',    'Capture perilune altitude (km)',  options.arrivalAltitude);
fprintf('  %-34s  %5.0f\n',    'Capture apolune altitude (km)',   options.arrivalApogeeAltitude);
fprintf('  %-34s  %5.1f deg\n','Capture inclination',             options.arrivalInclination);

fprintf('\nMission Budget\n');
fprintf('  %-34s  %9s  %9s\n', '',               'Direct', 'Bi-elliptic');
fprintf('  %-34s  %9s  %9s\n', '',               '------', '-----------');
fprintf('  %-34s  %7.3f    %7.3f  km/s\n', 'TLI delta-V', ...
    result_direct.details.dvTLI, result_bi.details.dvTLI);
fprintf('  %-34s  %7.3f    %7.3f  km/s\n', 'Total arrival delta-V', ...
    result_direct.details.dvCapture, result_bi.details.dvCapture);
fprintf('    %-32s  %7s    %7.3f  km/s\n', 'Burn 1 (equatorial LOI)', '--', result_bi.details.dvLOI);
fprintf('    %-32s  %7s    %7.3f  km/s\n', 'Burn 2 (plane change)', '--', result_bi.details.dvPlaneChange);
fprintf('    %-32s  %7s    %7.3f  km/s\n', 'Burn 3 (apoapsis trim)', '--', result_bi.details.dvApoapsisTrim);
fprintf('  %-34s  %7.0f    %7.0f  m/s\n',  'TCM budget', ...
    dv_tcm_dir*1000, dv_tcm_bi*1000);
fprintf('  %-34s  %7.3f    %7.3f  km/s\n', 'TOTAL delta-V (+TCM)', ...
    dv_total_dir, dv_total_bi);
fprintf('  %-34s  %7.2f    %7.2f  days\n', 'Time of flight', ...
    result_direct.tof / 86400, result_bi.tof / 86400);
fprintf('  %-34s  %7s    %7.0f  km\n', 'Intermediate apoapsis altitude', ...
    '--', result_bi.details.rBiEllipticApoapsis - Moon.radius);
fprintf('\n  %-34s  %+.3f km/s\n', 'Bi-elliptic savings (+TCM)', ...
    dv_total_dir - dv_total_bi);

fprintf('\nCapture Orbit Parameters\n');
fprintf('  Perilune altitude:           %5.0f km\n',    options.arrivalAltitude);
fprintf('  Apolune altitude:            %5.0f km\n',    options.arrivalApogeeAltitude);
fprintf('  Argument of periapsis:       %5.0f deg\n',   options.arrivalArgOfPeriapsis);
fprintf('  Semi-major axis:             %5.0f km\n',    a_cap);
fprintf('  Eccentricity:                %7.4f\n',        e_cap);
fprintf('  Orbital period:              %7.2f hours\n', T_cap_hr);
fprintf('  Speed at perilune:           %7.3f km/s\n',  v_peri);
fprintf('  Speed at apolune:            %7.3f km/s\n',  v_apo);

%% Trajectory plots — direct transfer (overview + body-centric)
plotPatchedConic(result_direct, Earth, Moon);
saveas(gcf, fullfile(outDir, 'lunar_south_pole_direct_transfer.png'));

%% Trajectory plots — bi-elliptic transfer (body-centric shows intermediate orbits)
plotPatchedConic(result_bi, Earth, Moon);
saveas(gcf, fullfile(outDir, 'lunar_south_pole_bielliptic_transfer.png'));

%% Delta-V comparison figure
figDV = figure('Name', 'Delta-V Comparison: Direct vs Bi-elliptic', 'NumberTitle', 'off', ...
       'Position', [100 100 620 420]);

labels  = {'TLI', 'LOI (equatorial)', 'Plane change', 'Apoapsis trim', 'TCM'};
dv_dir  = [result_direct.details.dvTLI,  result_direct.details.dvCapture, 0, 0, dv_tcm_dir];
dv_bi   = [result_bi.details.dvTLI, result_bi.details.dvLOI, ...
           result_bi.details.dvPlaneChange, result_bi.details.dvApoapsisTrim, dv_tcm_bi];

x = 1:5;
b1 = bar(x - 0.2, dv_dir, 0.35, 'FaceColor', [0.2 0.5 1.0]);
hold on;
b2 = bar(x + 0.2, dv_bi,  0.35, 'FaceColor', [0.9 0.45 0.1]);
hold off;

set(gca, 'XTick', x, 'XTickLabel', labels);
ylabel('\DeltaV (km/s)');
title('Burn Budget: Direct vs Bi-elliptic (90° polar capture)');
legend([b1 b2], ...
    sprintf('Direct  (total %.3f km/s)',      dv_total_dir), ...
    sprintf('Bi-elliptic  (total %.3f km/s)', dv_total_bi), ...
    'Location', 'northwest');
grid on;
text(0.98, 0.97, sprintf('Savings: %.3f km/s', dv_total_dir - dv_total_bi), ...
     'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
     'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.1 0.6 0.1]);
saveas(figDV, fullfile(outDir, 'lunar_south_pole_dv_comparison.png'));
