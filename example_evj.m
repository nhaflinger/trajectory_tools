% example_evj.m
% Earth -> Venus -> Jupiter gravity-assist trajectory (EVJ).
%
% Uses a Venus gravity assist to reduce departure C3 for a Jupiter mission.
% Searches the 2026-2040 window for the best opportunity and compares
% to a direct Earth -> Jupiter transfer.
%
% Jupiter is captured into a high-apojove elliptical orbit (~Juno-like,
% apojove ~8.1 million km) to keep the JOI delta-V small.

bodies = constants();
bSeq   = {bodies.Earth, bodies.Venus, bodies.Jupiter};

%% ---- Search for best EVJ launch window (2026-2040) ----
jdStart = julianDate(2026, 1, 1);
jdEnd   = julianDate(2040, 12, 31);

% Leg 1: Earth -> Venus  (~100-250 days)
% Leg 2: Venus -> Jupiter (~600-1200 days)
tofRanges = [100, 250; ...
             600, 1200];

searchOpts.nDepDates      = 60;
searchOpts.nTofPoints     = 30;
searchOpts.flybyAltitudes = 500;   % km above Venus surface

fprintf('Searching EVJ launch window 2026-2040...\n');
best = findBestFlybyWindow(bSeq, jdStart, jdEnd, tofRanges, searchOpts);

fprintf('  Best departure: %s\n', ...
    datestr(best.departureJD - 1721058.5, 'yyyy-mmm-dd'));
fprintf('  TOF leg 1 (E->V): %.0f days\n', best.tofDays(1));
fprintf('  TOF leg 2 (V->J): %.0f days\n', best.tofDays(2));
fprintf('  DV proxy:         %.3f km/s\n\n', best.deltaVProxy);

%% ---- Full trajectory at the best date ----
seqOpts = struct();
seqOpts.departureAltitude     = 200;       % km, LEO parking orbit
seqOpts.arrivalAltitude       = 500;       % km, Jupiter perijove altitude
seqOpts.arrivalApogeeAltitude = 8100000;   % km, highly elliptical capture
                                            % (~Juno-like, apojove ~116 R_J)
seqOpts.departureInclination  = 28.5;      % deg, KSC
seqOpts.arrivalInclination    = 0;         % deg
seqOpts.flybyAltitudes        = 500;       % km above Venus surface
seqOpts.atmosphereAltitudes   = 250;       % km, Venus atmosphere top

result = flybySequence(bSeq, best.departureJD, best.tofDays, seqOpts);

%% ---- Console output ----
dep_str = datestr(best.departureJD          - 1721058.5, 'yyyy-mmm-dd');
vfb_str = datestr(result.legs(1).arriveJD   - 1721058.5, 'yyyy-mmm-dd');
arr_str = datestr(result.legs(2).arriveJD   - 1721058.5, 'yyyy-mmm-dd');

fprintf('Earth -> Venus -> Jupiter  Gravity Assist\n');
fprintf('==========================================\n\n');

fprintf('Departure - Earth\n');
fprintf('  Date:                        %s\n',          dep_str);
fprintf('  Parking orbit altitude:      %5.0f km\n',    seqOpts.departureAltitude);
fprintf('  Parking orbit inclination:   %5.1f deg\n',   seqOpts.departureInclination);
fprintf('  C3:                          %6.2f km^2/s^2\n', result.details.C3);
fprintf('  v_inf departure:             %6.3f km/s\n',  result.details.vInfDepart);
fprintf('  Departure DV:                %6.3f km/s\n',  result.details.dvDeparture);

fprintf('\nLeg 1 - Earth to Venus\n');
fprintf('  Departure:                   %s\n',          dep_str);
fprintf('  Arrival:                     %s\n',          vfb_str);
fprintf('  Time of flight:              %5.0f days\n',  result.legs(1).tofDays);
fprintf('  v_inf at Venus (arriving):   %6.3f km/s\n',  result.legs(1).v_inf_arrive);

fprintf('\nVenus Flyby\n');
fb = result.flybys(1);
fprintf('  Date:                        %s\n',          vfb_str);
fprintf('  v_inf in:                    %6.3f km/s\n',  fb.v_inf_in);
fprintf('  v_inf out:                   %6.3f km/s\n',  fb.v_inf_out);
fprintf('  Deflection angle:            %6.2f deg\n',   fb.deflection);
fprintf('  Max achievable deflection:   %6.2f deg\n',   fb.maxDeflection);
fprintf('  Required periapsis altitude: %6.0f km\n',    fb.altitude);
fprintf('  Flyby feasible:              %s\n',          yesno(fb.isFeasible));
if fb.isPowered
    fprintf('  Powered flyby DV:            %6.3f km/s\n', fb.dvPowered);
else
    fprintf('  Type:                        Free gravity assist\n');
end

fprintf('\nLeg 2 - Venus to Jupiter\n');
fprintf('  Departure:                   %s\n',          vfb_str);
fprintf('  Arrival:                     %s\n',          arr_str);
fprintf('  Time of flight:              %5.0f days\n',  result.legs(2).tofDays);
fprintf('  v_inf at Jupiter (arriving): %6.3f km/s\n',  result.legs(2).v_inf_arrive);

fprintf('\nArrival - Jupiter\n');
fprintf('  Date:                        %s\n',          arr_str);
fprintf('  Perijove altitude:           %5.0f km\n',    seqOpts.arrivalAltitude);
fprintf('  Apojove altitude:            %.3e km\n',     seqOpts.arrivalApogeeAltitude);
fprintf('  Capture DV (JOI):            %6.3f km/s\n',  result.details.dvArrival);

fprintf('\nMission Budget\n');
fprintf('  Departure DV:                %6.3f km/s\n',  result.details.dvDeparture);
if result.details.dvPoweredFlybys > 1e-4
    fprintf('  Powered flyby DV:            %6.3f km/s\n', result.details.dvPoweredFlybys);
end
fprintf('  Arrival DV (JOI):            %6.3f km/s\n',  result.details.dvArrival);
fprintf('  TCM budget:                  %6.0f m/s\n',   result.details.dvTCM*1000);
fprintf('  Total DV (+TCM):             %6.3f km/s\n',  result.deltaV);
fprintf('  Total TOF:                   %6.0f days  (%.1f years)\n', ...
    result.tof, result.tof/365.25);

%% ---- Comparison: direct Earth -> Jupiter ----
fprintf('\nComparison - Direct Earth -> Jupiter\n');
fprintf('  (best launch in same window, no gravity assist)\n');

bestDirect = findBestLaunchDate(bodies.Earth, bodies.Jupiter, jdStart, jdEnd, ...
    linspace(600, 1400, 40)');

optDirect = struct();
optDirect.departureAltitude     = seqOpts.departureAltitude;
optDirect.arrivalAltitude       = seqOpts.arrivalAltitude;
optDirect.arrivalApogeeAltitude = seqOpts.arrivalApogeeAltitude;
optDirect.departureInclination  = seqOpts.departureInclination;
optDirect.arrivalInclination    = seqOpts.arrivalInclination;
optDirect.departureJD           = bestDirect.departureJD;
optDirect.tofDays               = bestDirect.tofDays;

resultDirect = patchedConicTransfer(bodies.Earth, bodies.Jupiter, optDirect);
dep_str_dir  = datestr(bestDirect.departureJD - 1721058.5, 'yyyy-mmm-dd');

fprintf('  Best departure:              %s\n',         dep_str_dir);
fprintf('  TOF:                         %5.0f days\n', bestDirect.tofDays);
fprintf('  Total DV (+TCM):             %6.3f km/s\n', resultDirect.deltaV);
fprintf('\n  EVJ gravity assist saves:    %+.3f km/s\n', ...
    resultDirect.deltaV - result.deltaV);

%% ---- Plot ----
plotFlybySequence(result, bSeq);

function s = yesno(tf)
    if tf, s = 'Yes'; else, s = 'No (powered flyby needed)'; end
end
