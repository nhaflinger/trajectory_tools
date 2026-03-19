function result = flybySequence(bodies, departureJD, tofDays, options)
%FLYBYSEQUENCE  Multi-leg patched-conic trajectory with gravity assists.
%
%   result = flybySequence(bodies, departureJD, tofDays)
%   result = flybySequence(bodies, departureJD, tofDays, options)
%
%   Chains N−1 Lambert legs through N bodies (first = departure, last =
%   arrival, middle = gravity-assist flyby bodies) and analyses the flyby
%   geometry at each intermediate body.
%
%   Inputs:
%     bodies       - 1×N cell array of body structs from constants()
%     departureJD  - Julian Date of departure from the first body
%     tofDays      - 1×(N−1) time-of-flight per leg (days)
%     options
%       .departureAltitude     parking orbit altitude at departure body (km) [200]
%       .arrivalAltitude       capture orbit periapsis altitude (km) [400]
%       .arrivalApogeeAltitude capture orbit apoapsis altitude (km) [= arrivalAltitude]
%       .departureInclination  parking orbit inclination (deg) [0]
%       .arrivalInclination    capture orbit inclination (deg) [0]
%       .flybyAltitudes        1×(N−2) flyby periapsis altitudes (km) [300 each]
%       .atmosphereAltitudes   1×(N−2) atmosphere-top altitudes (km) [0 each]
%       .transferTypes         1×(N−1) cell array 'type1'|'type2' per leg ['type1']
%
%   Output (result struct):
%     .legs        1×(N−1) struct array — per-leg Lambert solution
%     .flybys      1×(N−2) struct array — per-flyby geometry / feasibility
%     .deltaV      total ΔV: departure + powered flybys + arrival + TCM (km/s)
%     .deltaVBurns total ΔV without TCM (km/s)
%     .tof         total time of flight (days)
%     .details     budget breakdown and departure/arrival parameters
%
%   Per-leg fields (.legs(i)):
%     .departJD .arriveJD .tofDays
%     .r1_vec .v1_body .v1_transfer   — at leg departure
%     .r2_vec .v2_body .v2_transfer   — at leg arrival
%     .v_inf_depart_vec .v_inf_depart — departure hyperbolic excess (km/s)
%     .v_inf_arrive_vec .v_inf_arrive — arrival hyperbolic excess (km/s)
%
%   Per-flyby fields (.flybys(j)) — see gravityAssist.m:
%     .body .jd .v_inf_in .v_inf_out .deflection
%     .altitude .r_periapsis .isFeasible .isPowered .dvPowered .maxDeflection

if nargin < 4, options = struct(); end

N     = numel(bodies);
nLegs = N - 1;
nFB   = N - 2;
tofDays = tofDays(:)';

% --- Defaults ---
if ~isfield(options, 'departureAltitude'),    options.departureAltitude    = 200; end
if ~isfield(options, 'arrivalAltitude'),      options.arrivalAltitude      = 400; end
if ~isfield(options, 'arrivalApogeeAltitude')
    options.arrivalApogeeAltitude = options.arrivalAltitude;
end
if ~isfield(options, 'departureInclination'), options.departureInclination = 0;   end
if ~isfield(options, 'arrivalInclination'),   options.arrivalInclination   = 0;   end
if ~isfield(options, 'flybyAltitudes')
    options.flybyAltitudes = 300 * ones(1, max(nFB, 1));
end
if ~isfield(options, 'atmosphereAltitudes')
    options.atmosphereAltitudes = zeros(1, max(nFB, 1));
end
if ~isfield(options, 'transferTypes')
    options.transferTypes = repmat({'type1'}, 1, nLegs);
end
if ischar(options.transferTypes)
    options.transferTypes = repmat({options.transferTypes}, 1, nLegs);
end

muSun = constants().Sun.mu;

% --- Solve Lambert problem for each leg ---
legs = struct('departJD', {}, 'arriveJD', {}, 'tofDays', {}, ...
    'r1_vec', {}, 'v1_body', {}, 'v1_transfer', {}, ...
    'r2_vec', {}, 'v2_body', {}, 'v2_transfer', {}, ...
    'v_inf_depart_vec', {}, 'v_inf_depart', {}, ...
    'v_inf_arrive_vec', {}, 'v_inf_arrive', {});

jdAccum = departureJD;
for i = 1:nLegs
    jd_dep = jdAccum;
    jd_arr = jd_dep + tofDays(i);
    tof    = tofDays(i) * 86400;
    jdAccum = jd_arr;

    [r1_vec, v1_body] = orbitalState(bodies{i},   jd_dep);
    [r2_vec, v2_body] = orbitalState(bodies{i+1}, jd_arr);

    isLW = strcmpi(options.transferTypes{i}, 'type2');
    [v1t, v2t] = lambertSolver(r1_vec, r2_vec, tof, muSun, isLW);

    legs(i).departJD         = jd_dep;
    legs(i).arriveJD         = jd_arr;
    legs(i).tofDays          = tofDays(i);
    legs(i).r1_vec           = r1_vec;
    legs(i).v1_body          = v1_body;
    legs(i).v1_transfer      = v1t;
    legs(i).r2_vec           = r2_vec;
    legs(i).v2_body          = v2_body;
    legs(i).v2_transfer      = v2t;
    legs(i).v_inf_depart_vec = v1t - v1_body;
    legs(i).v_inf_depart     = norm(v1t - v1_body);
    legs(i).v_inf_arrive_vec = v2t - v2_body;
    legs(i).v_inf_arrive     = norm(v2t - v2_body);
end

% --- Analyse each flyby ---
% Build struct array from first assignment (pre-init with struct([]) causes
% "dissimilar structures" error when fb carries fields).
dv_powered = 0;
for j = 1:nFB
    fb_opts.atmosphereAltitude = options.atmosphereAltitudes(j);
    fb_opts.flybyAltitude      = options.flybyAltitudes(j);

    fb = gravityAssist(bodies{j+1}, ...
        legs(j).v_inf_arrive_vec, ...    % arriving on leg j
        legs(j+1).v_inf_depart_vec, ...  % departing on leg j+1
        fb_opts);
    fb.body = bodies{j+1};
    fb.jd   = legs(j).arriveJD;

    if j == 1
        flybys = fb;          % first element establishes the field set
    else
        flybys(j) = fb;       % subsequent elements extend the array
    end
    dv_powered = dv_powered + fb.dvPowered;
end
if nFB == 0
    flybys = struct([]);      % no flyby bodies — return empty struct array
end

% --- Departure burn (hyperbolic escape from parking orbit) ---
dep       = bodies{1};
r_park    = dep.radius + options.departureAltitude;
v_park    = sqrt(dep.mu / r_park);
v_inf_dep = legs(1).v_inf_depart;
v_hyp_dep = sqrt(v_inf_dep^2 + 2*dep.mu/r_park);
dv_dep    = combinedBurnDV(v_park, v_hyp_dep, options.departureInclination);
C3        = v_inf_dep^2;

% --- Arrival burn ---
arr       = bodies{N};
r_cap     = arr.radius + options.arrivalAltitude;
v_inf_arr = legs(nLegs).v_inf_arrive;
v_hyp_arr = sqrt(v_inf_arr^2 + 2*arr.mu/r_cap);
r_apo     = arr.radius + options.arrivalApogeeAltitude;
if r_apo > r_cap
    a_cap = (r_cap + r_apo) / 2;
    v_cap = sqrt(2*arr.mu/r_cap - arr.mu/a_cap);
else
    v_cap = sqrt(arr.mu / r_cap);
end
dv_arr = combinedBurnDV(v_hyp_arr, v_cap, options.arrivalInclination);

% --- TCM reserve ---
dv_tcm = max(0.010, min(0.050, 0.02 * v_inf_dep));

% --- Pack result ---
result.legs        = legs;
result.flybys      = flybys;
result.deltaVBurns = dv_dep + dv_powered + dv_arr;
result.deltaV      = dv_dep + dv_powered + dv_arr + dv_tcm;
result.tof         = sum(tofDays);
result.details     = struct( ...
    'dvDeparture',     dv_dep,    ...
    'dvPoweredFlybys', dv_powered, ...
    'dvArrival',       dv_arr,    ...
    'dvTCM',           dv_tcm,    ...
    'C3',              C3,        ...
    'vInfDepart',      v_inf_dep, ...
    'vInfArrive',      v_inf_arr, ...
    'departureJD',     departureJD, ...
    'tofDays',         tofDays,   ...
    'departureAltitude', options.departureAltitude, ...
    'arrivalAltitude',   options.arrivalAltitude);
end

% -------------------------------------------------------------------------
function dv = combinedBurnDV(v1, v2, inc_deg)
    dv = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(deg2rad(inc_deg)));
end
