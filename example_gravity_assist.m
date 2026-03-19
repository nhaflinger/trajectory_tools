% example_gravity_assist.m
% Earth -> Venus -> Jupiter gravity-assist trajectory.
%
% Demonstrates the patched-conic gravity-assist toolkit:
%   1. Grid search for the best launch window (2026-2040)
%   2. Full trajectory analysis at the optimal date
%   3. Comparison with a direct Earth -> Jupiter transfer
%   4. 3D heliocentric visualization

bodies = constants();
bSeq = {bodies.Earth, bodies.Venus, bodies.Jupiter};

%% ---- Search for best EVJ launch window (2026–2040) ----
jdStart = julianDate(2026, 1, 1);
jdEnd   = julianDate(2040, 12, 31);

% Leg 1: Earth -> Venus  (~100-250 days)
% Leg 2: Venus -> Jupiter (~600-1200 days)
tofRanges = [100, 250; ...
             600, 1200];

searchOpts.nDepDates  = 60;
searchOpts.nTofPoints = 30;
searchOpts.flybyAltitudes = 500;   % km above Venus surface

fprintf('Searching EVJ launch window 2026–2040...\n');
best = findBestFlybyWindow(bSeq, jdStart, jdEnd, tofRanges, searchOpts);

fprintf('  Best departure: %s\n', ...
    datestr(best.departureJD - 1721058.5, 'yyyy-mmm-dd'));
fprintf('  TOF leg 1 (E->V): %.0f days\n', best.tofDays(1));
fprintf('  TOF leg 2 (V->J): %.0f days\n', best.tofDays(2));
fprintf('  ΔV proxy:         %.3f km/s\n\n', best.deltaVProxy);

%% ---- Full trajectory at the best date ----
seqOpts = struct();
seqOpts.departureAltitude     = 200;       % km, LEO
seqOpts.arrivalAltitude       = 500;       % km, Jupiter perijove altitude
seqOpts.arrivalApogeeAltitude = 8100000;   % km, highly elliptical capture
                                            % (~Juno-like 53-day orbit, apojove ~116 R_J)
                                            % Captures with small ΔV; circularisation
                                            % done later via subsequent perijove burns.
seqOpts.departureInclination  = 28.5;      % deg, KSC
seqOpts.arrivalInclination    = 0;         % deg
seqOpts.flybyAltitudes        = 500;       % km above Venus surface
seqOpts.atmosphereAltitudes   = 250;       % km, Venus atmosphere top

result = flybySequence(bSeq, best.departureJD, best.tofDays, seqOpts);

%% ---- Console output ----
dep_str = datestr(best.departureJD - 1721058.5, 'yyyy-mmm-dd');
vfb_str = datestr(result.legs(1).arriveJD - 1721058.5, 'yyyy-mmm-dd');
arr_str = datestr(result.legs(2).arriveJD - 1721058.5, 'yyyy-mmm-dd');

fprintf('Earth -> Venus -> Jupiter  Gravity Assist\n');
fprintf('==========================================\n\n');

fprintf('Departure — Earth\n');
fprintf('  Date:                        %s\n', dep_str);
fprintf('  Parking orbit altitude:      %5.0f km\n',   seqOpts.departureAltitude);
fprintf('  Parking orbit inclination:   %5.1f deg\n',  seqOpts.departureInclination);
fprintf('  C3:                          %6.2f km²/s²\n', result.details.C3);
fprintf('  v∞ departure:                %6.3f km/s\n',   result.details.vInfDepart);
fprintf('  Departure ΔV:                %6.3f km/s\n',   result.details.dvDeparture);

fprintf('\nLeg 1 — Earth to Venus\n');
fprintf('  Departure:                   %s\n', dep_str);
fprintf('  Arrival:                     %s\n', vfb_str);
fprintf('  Time of flight:              %5.0f days\n',  result.legs(1).tofDays);
fprintf('  v∞ at Venus (arriving):      %6.3f km/s\n',  result.legs(1).v_inf_arrive);

fprintf('\nVenus Flyby\n');
fb = result.flybys(1);
fprintf('  Date:                        %s\n', vfb_str);
fprintf('  v∞ in:                       %6.3f km/s\n',  fb.v_inf_in);
fprintf('  v∞ out:                      %6.3f km/s\n',  fb.v_inf_out);
fprintf('  Deflection angle:            %6.2f deg\n',   fb.deflection);
fprintf('  Max achievable deflection:   %6.2f deg\n',   fb.maxDeflection);
fprintf('  Required periapsis altitude: %6.0f km\n',    fb.altitude);
fprintf('  Flyby feasible:              %s\n',           yesno(fb.isFeasible));
if fb.isPowered
    fprintf('  Powered flyby ΔV:            %6.3f km/s\n', fb.dvPowered);
else
    fprintf('  Type:                        Free gravity assist\n');
end

fprintf('\nLeg 2 — Venus to Jupiter\n');
fprintf('  Departure:                   %s\n', vfb_str);
fprintf('  Arrival:                     %s\n', arr_str);
fprintf('  Time of flight:              %5.0f days\n',  result.legs(2).tofDays);
fprintf('  v∞ at Jupiter (arriving):    %6.3f km/s\n',  result.legs(2).v_inf_arrive);

fprintf('\nArrival — Jupiter\n');
fprintf('  Date:                        %s\n', arr_str);
fprintf('  Capture periapsis altitude:  %5.0f km\n',   seqOpts.arrivalAltitude);
fprintf('  Capture ΔV:                  %6.3f km/s\n', result.details.dvArrival);

fprintf('\nMission Budget\n');
fprintf('  Departure ΔV:                %6.3f km/s\n',  result.details.dvDeparture);
if result.details.dvPoweredFlybys > 1e-4
    fprintf('  Powered flyby ΔV:            %6.3f km/s\n', result.details.dvPoweredFlybys);
end
fprintf('  Arrival ΔV:                  %6.3f km/s\n',  result.details.dvArrival);
fprintf('  TCM budget:                  %6.0f m/s\n',   result.details.dvTCM*1000);
fprintf('  Total ΔV (+TCM):             %6.3f km/s\n',  result.deltaV);
fprintf('  Total TOF:                   %6.0f days  (%.1f years)\n', ...
    result.tof, result.tof/365.25);

%% ---- Comparison: direct Earth -> Jupiter ----
fprintf('\nComparison — Direct Earth -> Jupiter\n');
fprintf('  (best launch in same window, no gravity assist)\n');

directOpts = struct();
directOpts.nDepDates  = 60;
directOpts.nTofPoints = 30;
tofDirectRange = [600, 1400];

jdEndDirect = julianDate(2040, 12, 31);
bestDirect = findBestLaunchDate(bodies.Earth, bodies.Jupiter, jdStart, jdEndDirect, ...
    linspace(tofDirectRange(1), tofDirectRange(2), directOpts.nTofPoints)');

directSeqOpts = struct();
directSeqOpts.departureAltitude    = seqOpts.departureAltitude;
directSeqOpts.arrivalAltitude      = seqOpts.arrivalAltitude;
directSeqOpts.arrivalApogeeAltitude = seqOpts.arrivalApogeeAltitude;
directSeqOpts.departureInclination = seqOpts.departureInclination;
directSeqOpts.arrivalInclination   = seqOpts.arrivalInclination;
directSeqOpts.departureJD          = bestDirect.departureJD;
directSeqOpts.tofDays              = bestDirect.tofDays;

resultDirect = patchedConicTransfer(bodies.Earth, bodies.Jupiter, directSeqOpts);
dep_str_dir  = datestr(bestDirect.departureJD - 1721058.5, 'yyyy-mmm-dd');

fprintf('  Best departure:              %s\n', dep_str_dir);
fprintf('  TOF:                         %5.0f days\n',  bestDirect.tofDays);
fprintf('  Total ΔV (+TCM):             %6.3f km/s\n',  resultDirect.deltaV);
fprintf('\n  EVJ gravity assist saves:    %+.3f km/s\n', ...
    resultDirect.deltaV - result.deltaV);

%% ---- Plots ----
plotFlybySequence(result, bSeq);

% ---- Helper ----
function s = yesno(tf)
    if tf, s = 'Yes'; else, s = 'No (powered flyby needed)'; end
end
