% example_lunar_transfer.m
% Earth -> Moon patched-conic transfer.
% Reports parking orbit parameters, TLI injection speed, transfer orbit
% geometry, v∞ at Moon SOI, capture ΔV, and TCM budget.

bodies = constants();
Earth  = bodies.Earth;
Moon   = bodies.Moon;

outDir = fullfile(fileparts(mfilename('fullpath')), 'output_lunar_transfer');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

options = struct();
options.departureAltitude    = 200;   % km, LEO parking orbit
options.arrivalAltitude      = 100;   % km, circular lunar orbit altitude
options.arrivalApogeeAltitude= 100;   % km, circular capture (= arrivalAltitude)
options.departureInclination = 28.5;  % deg, KSC parking orbit vs transfer plane
options.arrivalInclination   = 0;     % deg, equatorial lunar capture

result = patchedConicTransfer(Earth, Moon, options);

% ---- Derived transfer orbit quantities ----
r0           = result.details.r0;
a_xfr        = result.details.transferSemiMajor;
e_xfr        = 1 - r0 / a_xfr;                      % perigee = r0, apogee = Moon.a
T_xfr_days   = 2*pi * sqrt(a_xfr^3 / Earth.mu) / 86400;
v_park       = sqrt(Earth.mu / r0);                  % circular parking orbit speed
v_perigee    = sqrt(2*Earth.mu/r0 - Earth.mu/a_xfr); % TLI injection speed
v_inf_moon   = result.details.vInf;                  % hyperbolic excess at Moon SOI
r_cap        = result.details.rParkArrive;
v_perilune   = sqrt(v_inf_moon^2 + 2*Moon.mu/r_cap); % approach speed at perilune

% TCM budget: 2% of TLI ΔV, minimum 10 m/s
dv_tcm   = max(0.010, 0.02 * result.details.dvTLI);
dv_total = result.deltaV + dv_tcm;

% ---- Console output ----
fprintf('\nLunar Transfer  (Earth -> Moon)\n');
fprintf('================================\n\n');

fprintf('Departure — Earth\n');
fprintf('  Parking orbit altitude:    %5.0f km\n',    options.departureAltitude);
fprintf('  Parking orbit inclination: %5.1f deg\n',   options.departureInclination);
fprintf('  Parking orbit speed:       %7.3f km/s\n',  v_park);
fprintf('  TLI injection speed:       %7.3f km/s\n',  v_perigee);
fprintf('  TLI delta-V:               %7.3f km/s\n',  result.details.dvTLI);

fprintf('\nTransfer Orbit\n');
fprintf('  Semi-major axis:           %7.0f km\n',    a_xfr);
fprintf('  Eccentricity:              %7.4f\n',        e_xfr);
fprintf('  Full period:               %7.2f days\n',  T_xfr_days);
fprintf('  Time of flight (half):     %7.2f days\n',  result.tof / 86400);

fprintf('\nArrival — Moon\n');
fprintf('  v∞ at Moon SOI:            %7.3f km/s\n',  v_inf_moon);
fprintf('  Approach speed at perilune:%7.3f km/s\n',  v_perilune);
fprintf('  Capture orbit altitude:    %5.0f km\n',    options.arrivalAltitude);
fprintf('  Capture orbit inclination: %5.1f deg\n',   options.arrivalInclination);
fprintf('  Capture delta-V:           %7.3f km/s\n',  result.details.dvCapture);

fprintf('\nMission Budget\n');
fprintf('  TLI delta-V:               %7.3f km/s\n',  result.details.dvTLI);
fprintf('  Capture delta-V:           %7.3f km/s\n',  result.details.dvCapture);
fprintf('  TCM budget:                %7.0f m/s\n',   dv_tcm * 1000);
fprintf('  Total delta-V (+TCM):      %7.3f km/s\n',  dv_total);

% ---- Trajectory plots ----
plotPatchedConic(result, Earth, Moon);
saveas(gcf, fullfile(outDir, 'lunar_transfer.png'));
