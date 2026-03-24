% example_low_thrust.m
% Low-thrust orbit transfer examples — Phase 1.
%
% Demonstrates lowThrustSpiral() and plotLowThrustSpiral() for Earth orbit
% transfers and compares continuous low-thrust against impulsive burns.
%
% Cases:
%   1. LEO -> GEO  (no plane change)
%      Hall thruster: Isp = 1800 s, F = 300 mN
%   2. LEO -> GEO  (28.5 deg -> 0 deg plane change, KSC launch)
%      Edelbaum handles the combined manoeuvre optimally
%   3. LEO -> Lunar distance  (RK4 spiral + trajectory plot)
%      Ion engine: Isp = 3000 s, F = 100 mN
%   4. Propulsion trade study
%      Fixed mission (LEO -> GEO, no plane change): vary Isp and thrust
%      to compare propellant fraction and transfer time

bodies = constants();
Earth  = bodies.Earth;

% Common orbit radii
r_LEO  = Earth.radius + 200;      % km  200 km circular
r_GEO  = Earth.radius + 35786;    % km  GEO
r_Moon = bodies.Moon.a;            % km  ~384 400 km (lunar distance from Earth)

g0 = 9.80665e-3;  % km/s^2

outDir = fullfile(fileparts(mfilename('fullpath')), 'output_low_thrust');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('\n=================================================================\n');
fprintf('  Low-Thrust Orbit Transfer Examples\n');
fprintf('=================================================================\n\n');

%% ---- Impulsive Hohmann reference (LEO -> GEO, no plane change) ------
mu     = Earth.mu;
a_t    = (r_LEO + r_GEO) / 2;
v_LEO  = sqrt(mu / r_LEO);
v_GEO  = sqrt(mu / r_GEO);
v_peri = sqrt(2*mu/r_LEO  - mu/a_t);
v_apo  = sqrt(2*mu/r_GEO  - mu/a_t);
dv_hoh = (v_peri - v_LEO) + (v_GEO - v_apo);

fprintf('--- Impulsive Hohmann: LEO -> GEO (no plane change) ---\n');
fprintf('  Burn 1 (LEO):  %.4f km/s\n', v_peri - v_LEO);
fprintf('  Burn 2 (GEO):  %.4f km/s\n', v_GEO - v_apo);
fprintf('  Total:         %.4f km/s\n\n', dv_hoh);

%% ---- Case 1: Hall thruster, LEO -> GEO, no plane change --------------
fprintf('--- Case 1: Hall thruster (Isp=1800s, F=300mN), LEO->GEO ---\n');
opts1              = struct();
opts1.thrustN      = 0.300;   % N
opts1.isp          = 1800;    % s
opts1.wetMass      = 500;     % kg

res1 = lowThrustSpiral(Earth, r_LEO, r_GEO, opts1);
printSpiralResult(res1, dv_hoh);

%% ---- Case 2: Same thruster, add 28.5 deg plane change ---------------
fprintf('--- Case 2: Hall thruster + 28.5 deg plane change ---\n');
opts2 = opts1;
opts2.inclinationChangeDeg = 28.5;

res2 = lowThrustSpiral(Earth, r_LEO, r_GEO, opts2);

% Impulsive combined for comparison:
dv_imp_combined = sqrt(v_peri^2 + v_LEO^2 - 2*v_peri*v_LEO*cos(deg2rad(28.5))) ...
                + (v_GEO - v_apo);   % plane change at departure, circularise at GEO

fprintf('  Edelbaum dV:         %.4f km/s\n', res2.deltaV);
fprintf('  Impulsive combined:  %.4f km/s   (plane change at departure)\n', dv_imp_combined);
fprintf('  Edelbaum saving:     %.4f km/s   (distributed plane change benefit)\n', ...
    dv_imp_combined - res2.deltaV);
fprintf('  TOF:                 %.0f days\n', res2.tofDays);
fprintf('  Propellant:          %.1f kg  (%.1f%% of wet mass)\n', ...
    res2.propellantMass, 100*res2.propellantMass/opts2.wetMass);
fprintf('\n');

%% ---- Case 3: Ion engine, LEO -> Lunar distance ----------------------
fprintf('--- Case 3: Ion engine (Isp=3000s, F=100mN), LEO->Lunar distance ---\n');
opts3          = struct();
opts3.thrustN  = 0.100;   % N
opts3.isp      = 3000;    % s
opts3.wetMass  = 300;     % kg

res3 = lowThrustSpiral(Earth, r_LEO, r_Moon, opts3);
printSpiralResult(res3, NaN);

% Plot the spiral
fig3 = plotLowThrustSpiral(res3, Earth);
saveas(fig3, fullfile(outDir, 'low_thrust_spiral_leo_lunar.png'));

%% ---- Case 4: Propulsion trade study ---------------------------------
fprintf('--- Case 4: Propulsion trade — LEO->GEO, no plane change ---\n');
fprintf('  %-20s %8s %8s %10s %12s %10s\n', ...
    'Propulsion type', 'Isp (s)', 'F (mN)', 'dV (km/s)', 'Prop (kg)', 'TOF (days)');
fprintf('  %s\n', repmat('-', 1, 72));

tradeIsp    = [450,   1800,  3000,  5000];
tradeThrust = [5000,   300,   100,    50];   % mN, scaled to be roughly comparable
tradeNames  = {'Chemical (biprop)', 'Hall thruster', 'Ion engine (Xe)', 'Ion engine (high-Isp)'};

for k = 1:numel(tradeIsp)
    ot         = struct();
    ot.thrustN = tradeThrust(k) * 1e-3;   % mN -> N
    ot.isp     = tradeIsp(k);
    ot.wetMass = 500;
    rt = lowThrustSpiral(Earth, r_LEO, r_GEO, ot);
    fprintf('  %-20s %8.0f %8.0f %10.4f %12.1f %10.1f\n', ...
        tradeNames{k}, tradeIsp(k), tradeThrust(k), ...
        rt.deltaV, rt.propellantMass, rt.tofDays);
end
fprintf('\n');
fprintf('  Note: impulsive Hohmann dV = %.4f km/s (Oberth benefit, higher dV for spiral)\n\n', dv_hoh);

%% ---- Plot Case 1 (GEO spiral) ---------------------------------------
fig1 = plotLowThrustSpiral(res1, Earth);
saveas(fig1, fullfile(outDir, 'low_thrust_spiral_leo_geo.png'));



%% ---- Helper ---------------------------------------------------------
function printSpiralResult(res, dv_impulsive)
fprintf('  Edelbaum dV:  %.4f km/s', res.deltaV);
if ~isnan(dv_impulsive)
    fprintf('   (impulsive: %.4f km/s  =>  Oberth deficit: +%.4f km/s)', ...
        dv_impulsive, res.deltaV - dv_impulsive);
end
fprintf('\n');
fprintf('  TOF:          %.0f days\n', res.tofDays);
fprintf('  Propellant:   %.1f kg  (%.1f%% of %.0f kg wet mass)\n', ...
    res.propellantMass, 100*res.propellantMass/res.details.wetMass, res.details.wetMass);
fprintf('  Final mass:   %.1f kg\n\n', res.finalMass);
end
