% example_emej_europa_clipper.m
% Europa Clipper-style Earth -> Mars -> Earth -> Jupiter trajectory (EMEJ).
%
% Models the MEGA (Mars-Earth Gravity Assist) architecture used by NASA's
% Europa Clipper mission (launched October 2024).
%
% Actual Europa Clipper reference timeline:
%   Launch:       14 Oct 2024
%   Mars flyby:   01 Mar 2025  (~138 days after launch)
%   Earth flyby:  03 Dec 2026  (~641 days after Mars flyby)
%   Jupiter arr:  11 Apr 2030  (~1226 days after Earth flyby, ~2005 days total)
%   Actual C3:    ~6.2 km^2/s^2  (departure DV ~3.5 km/s)
%   Actual JOI:   ~0.9 km/s
%   Actual total: ~4.5 km/s
%
% Transfer type assignment — the sign of dm in the BM&W Lambert formula
% determines the arc direction based on the CCW angle between the two
% position vectors:
%   CCW angle < 180 deg  ->  type1 = prograde short arc
%   CCW angle > 180 deg  ->  type2 = prograde long arc  (type1 = retrograde)
%
% For the MEGA M->E leg the CCW angle from Mars (Mar 2025, ~138 deg
% ecliptic lon) to Earth (Dec 2026, ~70 deg) is about 292 deg > 180 deg,
% so type2 gives the correct prograde arc that sweeps ~292 deg around the
% sun before re-encountering Earth.
%
% Leg types: E->M = type1,  M->E = type2,  E->J = type1
%
% Departure and Jupiter arrival dates are fixed to the actual mission
% dates.  Mars and Earth flyby epochs are optimised over a grid then
% refined with fminsearch to minimise total delta-V.

bodies = constants();
bSeq   = {bodies.Earth, bodies.Mars, bodies.Earth, bodies.Jupiter};

%% ---- Fixed boundary dates (actual EC mission) ----
jd_launch  = julianDate(2024, 10, 14);
jd_jup_arr = julianDate(2030,  4, 11);
t_total    = jd_jup_arr - jd_launch;   % days

dep_str = datestr(jd_launch   - 1721058.5, 'yyyy-mmm-dd');
arr_str = datestr(jd_jup_arr  - 1721058.5, 'yyyy-mmm-dd');

fprintf('Europa Clipper EMEJ Trajectory (MEGA architecture)\n');
fprintf('====================================================\n');
fprintf('Departure:       %s  (fixed — actual EC launch date)\n', dep_str);
fprintf('Jupiter arrival: %s  (fixed — actual EC arrival date)\n', arr_str);
fprintf('Total TOF:       %.0f days  (%.1f years)\n\n', t_total, t_total/365.25);

%% ---- Trajectory options ----
seqOpts = struct();
seqOpts.departureAltitude     = 200;
seqOpts.arrivalAltitude       = 500;
seqOpts.arrivalApogeeAltitude = 14000000;
seqOpts.departureInclination  = 28.5;
seqOpts.arrivalInclination    = 0;
seqOpts.flybyAltitudes        = [300, 300];
seqOpts.atmosphereAltitudes   = [50,  100];
seqOpts.transferTypes         = {'type1', 'type2', 'type1'};

%% ---- Grid search over flyby epochs ----
% Free variables: t_m = days after launch for Mars flyby,
%                 t_e = days after launch for Earth flyby.
% Leg 3 TOF = t_total - t_e (Jupiter arrival is fixed).
% Objective: minimise total delta-V.
%
% type2 for M->E gives the prograde >180-deg arc only when the CCW angle
% from Mars's position to Earth's position is >180 deg.  For other dates
% the arc becomes retrograde (high DV); the total-DV minimisation
% naturally avoids those regions.
fprintf('Grid search over flyby epochs (type1 / type2 / type1)...\n');

t_m_grid = linspace( 90, 200, 30);   % E->M TOF (days after launch)
t_e_grid = linspace(550, 950, 36);   % days after launch for Earth flyby

bestDV_grid   = Inf;
bestTofs_grid = [];

for im = 1:numel(t_m_grid)
    t_m = t_m_grid(im);
    for ie = 1:numel(t_e_grid)
        t_e = t_e_grid(ie);
        if t_e < t_m + 200, continue; end   % M->E needs time for long arc
        t3  = t_total - t_e;
        if t3 < 700, continue; end
        tofs = [t_m,  t_e - t_m,  t3];
        try
            res = flybySequence(bSeq, jd_launch, tofs, seqOpts);
            if res.deltaV < bestDV_grid
                bestDV_grid   = res.deltaV;
                bestTofs_grid = tofs;
            end
        catch
        end
    end
end

fprintf('  Best total DV (grid): %.3f km/s\n', bestDV_grid);
fprintf('  Mars flyby at day:  %.1f  (%s)\n', bestTofs_grid(1), ...
    datestr(jd_launch + bestTofs_grid(1) - 1721058.5, 'yyyy-mmm-dd'));
fprintf('  Earth flyby at day: %.1f  (%s)\n', sum(bestTofs_grid(1:2)), ...
    datestr(jd_launch + sum(bestTofs_grid(1:2)) - 1721058.5, 'yyyy-mmm-dd'));
fprintf('\n');

%% ---- fminsearch refinement ----
fprintf('Refining with fminsearch...\n');

x0   = [bestTofs_grid(1),  sum(bestTofs_grid(1:2))];
obj  = @(x) emej_obj(x, jd_launch, t_total, bSeq, seqOpts);
fopt = optimset('TolFun', 1e-10, 'TolX', 0.05, 'MaxFunEvals', 8000, 'Display', 'off');

[x_ref, dv_ref] = fminsearch(obj, x0, fopt);

t_m_ref  = x_ref(1);
t_e_ref  = x_ref(2);
tofs_ref = [t_m_ref,  t_e_ref - t_m_ref,  t_total - t_e_ref];
bestRes  = flybySequence(bSeq, jd_launch, tofs_ref, seqOpts);

mfb_str = datestr(jd_launch + t_m_ref - 1721058.5, 'yyyy-mmm-dd');
efb_str = datestr(jd_launch + t_e_ref - 1721058.5, 'yyyy-mmm-dd');

fprintf('  Best total DV (refined): %.3f km/s\n', dv_ref);
fprintf('  Mars flyby:  day %4.0f  (%s)\n', t_m_ref, mfb_str);
fprintf('  Earth flyby: day %4.0f  (%s)\n', t_e_ref, efb_str);
fprintf('  (Actual EC: Mars day 138 [2025-Mar-01], Earth day 779 [2026-Dec-03])\n\n');

%% ---- Full mission report ----
fprintf('Earth -> Mars -> Earth -> Jupiter  (EMEJ / MEGA)\n');
fprintf('==================================================\n\n');

fprintf('Departure - Earth  (%s)\n', dep_str);
fprintf('  C3:                          %6.2f km^2/s^2\n', bestRes.details.C3);
fprintf('  v_inf departure:             %6.3f km/s\n',     bestRes.details.vInfDepart);
fprintf('  Departure DV:                %6.3f km/s\n',     bestRes.details.dvDeparture);
fprintf('  (Actual EC departure DV:     ~3.50 km/s)\n');

fprintf('\nLeg 1 - Earth to Mars  (%s -> %s,  %.0f days)  [type1, prograde]\n', ...
    dep_str, mfb_str, bestRes.legs(1).tofDays);

fprintf('\nMars Flyby  (%s)\n', mfb_str);
printFlyby(bestRes.flybys(1));

fprintf('\nLeg 2 - Mars to Earth  (%s -> %s,  %.0f days)  [type2, prograde >180 deg]\n', ...
    mfb_str, efb_str, bestRes.legs(2).tofDays);

fprintf('\nEarth Flyby  (%s)\n', efb_str);
printFlyby(bestRes.flybys(2));

fprintf('\nLeg 3 - Earth to Jupiter  (%s -> %s,  %.0f days)  [type1, prograde]\n', ...
    efb_str, arr_str, bestRes.legs(3).tofDays);
fprintf('  v_inf at Jupiter (arriving): %6.3f km/s\n', bestRes.legs(3).v_inf_arrive);

fprintf('\nArrival - Jupiter  (%s)\n', arr_str);
fprintf('  Perijove altitude:           %5.0f km\n', seqOpts.arrivalAltitude);
fprintf('  Apojove altitude:            %.3e km  (~%.0f R_J)\n', ...
    seqOpts.arrivalApogeeAltitude, seqOpts.arrivalApogeeAltitude/bodies.Jupiter.radius);
fprintf('  Capture DV (JOI):            %6.3f km/s\n', bestRes.details.dvArrival);
fprintf('  (Actual EC JOI:              ~0.90 km/s)\n');

fprintf('\nMission Budget\n');
fprintf('  Departure DV:                %6.3f km/s\n',  bestRes.details.dvDeparture);
if bestRes.details.dvPoweredFlybys > 1e-4
    fprintf('  Powered flyby DV:            %6.3f km/s\n', bestRes.details.dvPoweredFlybys);
else
    fprintf('  Gravity assists:             Free (no powered flyby needed)\n');
end
fprintf('  Arrival DV (JOI):            %6.3f km/s\n',  bestRes.details.dvArrival);
fprintf('  TCM budget:                  %6.0f m/s\n',   bestRes.details.dvTCM*1000);
fprintf('  Total DV (+TCM):             %6.3f km/s\n',  bestRes.deltaV);
fprintf('  (Actual EC total:            ~4.50 km/s)\n');
fprintf('  Total TOF:                   %6.0f days  (%.1f years)\n', ...
    bestRes.tof, bestRes.tof/365.25);

%% ---- Comparison: direct Earth -> Jupiter, same departure date ----
fprintf('\nComparison - Direct Earth -> Jupiter  (same %s departure)\n', dep_str);

tofVec    = linspace(600, 1400, 50)';
bestDVdir = Inf;  bestTOFdir = NaN;
optDirect = struct('departureAltitude', seqOpts.departureAltitude, ...
    'arrivalAltitude', seqOpts.arrivalAltitude, ...
    'arrivalApogeeAltitude', seqOpts.arrivalApogeeAltitude, ...
    'departureInclination', seqOpts.departureInclination, ...
    'arrivalInclination', seqOpts.arrivalInclination, ...
    'departureJD', jd_launch);
for k = 1:numel(tofVec)
    optDirect.tofDays = tofVec(k);
    try
        rd = patchedConicTransfer(bodies.Earth, bodies.Jupiter, optDirect);
        if rd.deltaV < bestDVdir
            bestDVdir = rd.deltaV;  bestTOFdir = tofVec(k);
        end
    catch, end
end
optDirect.tofDays = bestTOFdir;
resultDirect = patchedConicTransfer(bodies.Earth, bodies.Jupiter, optDirect);
fprintf('  Best TOF:        %5.0f days  (arrival %s)\n', bestTOFdir, ...
    datestr(jd_launch + bestTOFdir - 1721058.5, 'yyyy-mmm-dd'));
fprintf('  Total DV:        %6.3f km/s\n', resultDirect.deltaV);
fprintf('\n  EMEJ saves vs direct EJ:  %+.3f km/s\n', ...
    resultDirect.deltaV - bestRes.deltaV);

%% ---- Plot ----
plotFlybySequence(bestRes, bSeq);

%% ---- Local helpers ----
function dv = emej_obj(x, jd_launch, t_total, bSeq, seqOpts)
    t_mars  = x(1);
    t_earth = x(2);
    t3      = t_total - t_earth;
    if t_mars < 60 || t_earth < t_mars + 150 || t3 < 600
        dv = 1e6; return
    end
    tofs = [t_mars,  t_earth - t_mars,  t3];
    try
        res = flybySequence(bSeq, jd_launch, tofs, seqOpts);
        dv  = res.deltaV;
    catch
        dv = 1e6;
    end
end

function printFlyby(fb)
    fprintf('  v_inf in:                    %6.3f km/s\n', fb.v_inf_in);
    fprintf('  v_inf out:                   %6.3f km/s\n', fb.v_inf_out);
    fprintf('  Deflection angle:            %6.2f deg\n',  fb.deflection);
    fprintf('  Max achievable deflection:   %6.2f deg\n',  fb.maxDeflection);
    fprintf('  Required periapsis altitude: %6.0f km\n',   fb.altitude);
    fprintf('  Flyby feasible:              %s\n',         yesno(fb.isFeasible));
    if fb.isPowered
        fprintf('  Powered flyby DV:            %6.3f km/s\n', fb.dvPowered);
    else
        fprintf('  Type:                        Free gravity assist\n');
    end
end

function s = yesno(tf)
    if tf, s = 'Yes'; else, s = 'No (powered flyby needed)'; end
end
