function result = patchedConicTransfer(departBody, arrivalBody, options)
%PATCHEDCONICTRANSFER Simple patched-conic trajectory approximations
%   result = patchedConicTransfer(departBody, arrivalBody, options)
%
%   Calculates a patched conic transfer between two solar system bodies.
%   This is a simplifed model assuming circular coplanar orbits and
%   uses Hohmann-like transfers for interplanetary legs. For lunar transfers
%   it uses a two-phase Earth-centered transfer (parking orbit -> trans-lunar
%   injection -> lunar capture).
%
%   Inputs:
%       departBody  - struct for departure body (from constants.m)
%       arrivalBody - struct for arrival body (from constants.m)
%       options     - struct with fields:
%                     * departureAltitude (km) [default: 200]
%                     * arrivalAltitude (km) [default: 200]
%                     * departureInclination (deg) inclination of parking orbit
%                         relative to the transfer plane [default: 0]
%                     * arrivalInclination (deg) inclination of capture orbit
%                         relative to the transfer plane [default: 0]
%                     * arrivalApogeeAltitude (km) apogee altitude of elliptical
%                         capture orbit; defaults to arrivalAltitude (circular)
%                     * nRevs (int) for Lambert solver [default: 0]
%
%   Outputs:
%       result - struct with fields:
%                * deltaV (km/s) total approximation
%                * tof (s) time-of-flight
%                * phaseAngle (deg) required phase angle at departure
%                * details (struct) with breakdowns

% Backwards-compatible input parsing (pre-R2019b)
if nargin < 3 || isempty(options)
    options = struct();
end

if ~isfield(options, 'departureAltitude') || isempty(options.departureAltitude)
    options.departureAltitude = 200;
end

if ~isfield(options, 'arrivalAltitude') || isempty(options.arrivalAltitude)
    options.arrivalAltitude = 200;
end

if ~isfield(options, 'nRevs') || isempty(options.nRevs)
    options.nRevs = 0;
end

if ~isfield(options, 'departureInclination') || isempty(options.departureInclination)
    options.departureInclination = 0;
end

if ~isfield(options, 'arrivalInclination') || isempty(options.arrivalInclination)
    options.arrivalInclination = 0;
end

if ~isfield(options, 'arrivalApogeeAltitude') || isempty(options.arrivalApogeeAltitude)
    options.arrivalApogeeAltitude = options.arrivalAltitude;
end

if ~isfield(options, 'arrivalArgOfPeriapsis') || isempty(options.arrivalArgOfPeriapsis)
    options.arrivalArgOfPeriapsis = 0;  % deg; 90 = perilune above north pole
end

if ~isfield(options, 'transferMode') || isempty(options.transferMode)
    options.transferMode = 'direct';
end

% Lambert options (interplanetary only)
if ~isfield(options, 'departureJD') || isempty(options.departureJD)
    options.departureJD = [];
end
if ~isfield(options, 'tofDays') || isempty(options.tofDays)
    options.tofDays = [];
end
if ~isfield(options, 'transferType') || isempty(options.transferType)
    options.transferType = 'type1';   % 'type1' = short-way, 'type2' = long-way
end

% biEllipticApoapsisAltitude: altitude (km) of intermediate equatorial orbit
% apoapsis for the plane-change burn.  0 = use Moon SOI (maximum savings).
if ~isfield(options, 'biEllipticApoapsisAltitude') || isempty(options.biEllipticApoapsisAltitude)
    options.biEllipticApoapsisAltitude = 0;
end

% Quick path for Earth->Moon lunar transfer
if strcmpi(departBody.name, 'Earth') && strcmpi(arrivalBody.name, 'Moon')
    if strcmpi(options.transferMode, 'biElliptic')
        result = lunarTransferBiElliptic(options.departureAltitude, options.arrivalAltitude, options);
    else
        result = lunarTransfer(options.departureAltitude, options.arrivalAltitude, options);
    end
    return;
end

% Otherwise treat it as interplanetary transfer around the Sun
result = interplanetaryTransfer(departBody, arrivalBody, options);
end

function res = lunarTransfer(h0, h1, options)
% LUNARTRANSFER Simple patched conic for Earth->Moon
bodies = constants();
Earth = bodies.Earth;
Moon = bodies.Moon;

% 1) Parking orbit around Earth
r0 = Earth.radius + h0; % km
v0 = sqrt(Earth.mu / r0);

% 2) Transfer ellipse to lunar distance (approx circular)
a_transfer = (r0 + Moon.a)/2;

% TLI burn: combined speed change + plane change from parking orbit to
% transfer ellipse. departureInclination is the angle between the parking
% orbit plane and the lunar transfer plane.
v_perigee = sqrt(2*Earth.mu/r0 - Earth.mu/a_transfer);
dv_tli = combinedBurnDV(v0, v_perigee, options.departureInclination);

% 3) Arrival at lunar distance: velocity in transfer orbit
v_arrival = sqrt(2*Earth.mu/Moon.a - Earth.mu/a_transfer);

% 4) Capture into lunar orbit at h1
r_peri_cap = Moon.radius + h1;
% Patched conic: v_inf relative to Moon = difference between Moon's circular
% velocity and spacecraft apogee velocity (both Earth-centered, co-aligned
% for a prograde transfer).
v_moon = sqrt(Earth.mu / Moon.a);
v_inf  = abs(v_moon - v_arrival);
% In Moon-centered frame, arrival speed at perilune.
v_perilune = sqrt(v_inf^2 + 2*Moon.mu/r_peri_cap);
% Target speed at perilune: elliptical orbit if arrivalApogeeAltitude > h1,
% circular otherwise.
r_apo_cap = Moon.radius + options.arrivalApogeeAltitude;
if r_apo_cap > r_peri_cap
    a_cap       = (r_peri_cap + r_apo_cap) / 2;
    v_cap_peri  = sqrt(2*Moon.mu/r_peri_cap - Moon.mu/a_cap);
else
    v_cap_peri  = sqrt(Moon.mu / r_peri_cap);
    a_cap       = r_peri_cap;
end
dv_capture = combinedBurnDV(v_perilune, v_cap_peri, options.arrivalInclination);

% Time of flight (half ellipse)
tof = pi*sqrt(a_transfer^3/Earth.mu);

res = struct();
res.deltaV = dv_tli + dv_capture;
res.tof = tof;
% Estimate required phase angle (radians) so Moon arrives near apogee
nMoon = sqrt(Earth.mu / Moon.a^3);
phaseAngle = rad2deg(pi - nMoon * tof);

res.phaseAngle = phaseAngle;
res.details = struct('dvTLI', dv_tli, ...
                         'dvCapture', dv_capture, ...
                         'vInf', v_inf, ...
                         'transferSemiMajor', a_transfer, ...
                         'r0', r0, ...
                         'rApogee', Moon.a, ...
                         'rParkArrive', r_peri_cap, ...
                         'rApoArrive',  r_apo_cap, ...
                         'aCaptureOrbit', a_cap, ...
                         'departureInclination', options.departureInclination, ...
                         'arrivalInclination',   options.arrivalInclination, ...
                         'arrivalArgOfPeriapsis', options.arrivalArgOfPeriapsis, ...
                         'tof', tof);
end

function res = interplanetaryTransfer(departBody, arrivalBody, options)
% INTERPLANETARYTRANSFER Patched-conic interplanetary transfer.
% Uses Lambert solver when departureJD + tofDays are provided;
% falls back to Hohmann approximation otherwise.

muSun = constants().Sun.mu;

% ---- Heliocentric leg ----
if ~isempty(options.departureJD) && ~isempty(options.tofDays)
    % Lambert-based transfer with true Keplerian planet positions
    jd_dep    = options.departureJD;
    jd_arr    = jd_dep + options.tofDays;
    tof       = options.tofDays * 86400;
    isLongWay = strcmpi(options.transferType, 'type2');

    [r1_vec, v1_body] = orbitalState(departBody, jd_dep);
    [r2_vec, v2_body] = orbitalState(arrivalBody, jd_arr);
    [v1t, v2t]        = lambertSolver(r1_vec, r2_vec, tof, muSun, isLongWay);

    v_inf_dep = norm(v1t - v1_body);
    v_inf_arr = norm(v2_body - v2t);
    r1        = norm(r1_vec);
    r2        = norm(r2_vec);
    a_t       = 1 / (2/r1 - dot(v1t, v1t)/muSun);
    phase_deg = 0;   % not a meaningful scalar for Lambert

    has_lambert = true;
else
    % Hohmann-style fallback (circular coplanar)
    r1 = departBody.a;   r2 = arrivalBody.a;
    a_t   = (r1 + r2) / 2;
    v_peri = sqrt(2*muSun/r1 - muSun/a_t);
    v_apo  = sqrt(2*muSun/r2 - muSun/a_t);
    v_inf_dep = abs(v_peri - sqrt(muSun/r1));
    v_inf_arr = abs(sqrt(muSun/r2) - v_apo);
    tof       = pi * sqrt(a_t^3 / muSun);
    phase_deg = rad2deg(pi*(1 - sqrt((2*r2)/(r1+r2))));

    r1_vec = [];  v1_body = [];  v1t = [];
    r2_vec = [];  v2_body = [];  v2t = [];
    has_lambert = false;
end

% ---- Body-centric burns (same for both paths) ----
r_park    = departBody.radius + options.departureAltitude;
v_park    = sqrt(departBody.mu / r_park);
v_hyp_dep = sqrt(v_inf_dep^2 + 2*departBody.mu/r_park);
dv_dep    = combinedBurnDV(v_park, v_hyp_dep, options.departureInclination);

r_cap     = arrivalBody.radius + options.arrivalAltitude;
r_apo_arr = arrivalBody.radius + options.arrivalApogeeAltitude;
if r_apo_arr > r_cap
    a_cap = (r_cap + r_apo_arr) / 2;
    v_cap = sqrt(2*arrivalBody.mu/r_cap - arrivalBody.mu/a_cap);
else
    v_cap = sqrt(arrivalBody.mu / r_cap);
end
v_hyp_arr = sqrt(v_inf_arr^2 + 2*arrivalBody.mu/r_cap);
dv_arr    = combinedBurnDV(v_hyp_arr, v_cap, options.arrivalInclination);

% ---- C3 and TCM budget (suggestions 4 & 5) ----
C3     = v_inf_dep^2;                              % km^2/s^2
dv_tcm = max(0.010, min(0.050, 0.02*v_inf_dep));   % 2% of v_inf_dep, 10-50 m/s

% ---- Pack result ----
res           = struct();
res.deltaV    = dv_dep + dv_arr + dv_tcm;
res.deltaVBurns = dv_dep + dv_arr;
res.tof       = tof;
res.phaseAngle = phase_deg;
res.details   = struct( ...
    'dvDeparture', dv_dep, ...
    'dvArrival',   dv_arr, ...
    'dvTCM',       dv_tcm, ...
    'C3',          C3, ...
    'vInfDepart',  v_inf_dep, ...
    'vInfArrive',  v_inf_arr, ...
    'tof',         tof, ...
    'r1',          r1, ...
    'r2',          r2, ...
    'aTransfer',   a_t, ...
    'rParkDepart', r_park, ...
    'rParkArrive', r_cap, ...
    'departureInclination', options.departureInclination, ...
    'arrivalInclination',   options.arrivalInclination);

if has_lambert
    res.details.r1_vec     = r1_vec;
    res.details.v1_transfer = v1t;
    res.details.r2_vec     = r2_vec;
    res.details.v2_transfer = v2t;
    res.details.v1_body    = v1_body;
    res.details.v2_body    = v2_body;
    res.details.departureJD = options.departureJD;
    res.details.tofDays    = options.tofDays;
end
end

function res = lunarTransferBiElliptic(h0, h1, options)
%LUNARTRANSFERBIELLIPTIC Bi-elliptic plane-change Earth->Moon transfer.
%
%   Arrival sequence (3 burns):
%     Burn 1 (perilune):  hyperbola -> equatorial ellipse, NO plane change
%     Burn 2 (apolune):   pure plane change where spacecraft is slowest
%     Burn 3 (perilune):  lower apoapsis from intermediate to final value
%
%   Key identity: dv_arr = (v_hyp_peri - v_final_peri) + 2*v_bi_apo*sin(i/2)
%   The v_bi_peri terms in burns 1 and 3 cancel exactly, so savings grow
%   monotonically with intermediate apoapsis.  Default uses Moon SOI.

bodies = constants();
Earth  = bodies.Earth;
Moon   = bodies.Moon;

% Moon SOI relative to Earth (used as default upper bound)
r_soi_moon = Moon.a * (Moon.mu / Earth.mu)^(2/5);  % ~66 183 km

% TLI — identical to direct case
r0         = Earth.radius + h0;
v0         = sqrt(Earth.mu / r0);
a_transfer = (r0 + Moon.a) / 2;
v_perigee  = sqrt(2*Earth.mu/r0 - Earth.mu/a_transfer);
dv_tli     = combinedBurnDV(v0, v_perigee, options.departureInclination);

% Arrival v_inf at Moon
v_arrival  = sqrt(2*Earth.mu/Moon.a - Earth.mu/a_transfer);
v_moon     = sqrt(Earth.mu / Moon.a);
v_inf      = abs(v_moon - v_arrival);

r_peri     = Moon.radius + h1;
v_hyp_peri = sqrt(v_inf^2 + 2*Moon.mu/r_peri);

r_apo_final = Moon.radius + options.arrivalApogeeAltitude;
a_final     = (r_peri + r_apo_final) / 2;
v_final_peri = sqrt(2*Moon.mu/r_peri - Moon.mu/a_final);

% Intermediate apoapsis: user-specified or Moon SOI (optimal)
if options.biEllipticApoapsisAltitude > 0
    r_bi_apo = Moon.radius + options.biEllipticApoapsisAltitude;
    r_bi_apo = max(r_bi_apo, r_apo_final);   % must be >= final orbit apoapsis
else
    r_bi_apo = r_soi_moon;
end

[dv_loi, dv_plane, dv_trim] = biEllipticCapture( ...
    Moon, v_hyp_peri, r_peri, r_bi_apo, r_apo_final, options.arrivalInclination);
dv_arr = dv_loi + dv_plane + dv_trim;

% Capture orbit for consistent reporting (same as direct)
if r_apo_final > r_peri
    a_cap = a_final;
else
    a_cap = r_peri;
end

tof        = pi * sqrt(a_transfer^3 / Earth.mu);
nMoon      = sqrt(Earth.mu / Moon.a^3);
phaseAngle = rad2deg(pi - nMoon * tof);

res           = struct();
res.deltaV    = dv_tli + dv_arr;
res.tof       = tof;
res.phaseAngle = phaseAngle;
res.details   = struct( ...
    'dvTLI',                  dv_tli, ...
    'dvCapture',              dv_arr, ...
    'dvLOI',                  dv_loi, ...
    'dvPlaneChange',          dv_plane, ...
    'dvApoapsisTrim',         dv_trim, ...
    'vInf',                   v_inf, ...
    'transferSemiMajor',      a_transfer, ...
    'r0',                     r0, ...
    'rApogee',                Moon.a, ...
    'rParkArrive',            r_peri, ...
    'rApoArrive',             r_apo_final, ...
    'rBiEllipticApoapsis',    r_bi_apo, ...
    'aCaptureOrbit',          a_cap, ...
    'departureInclination',   options.departureInclination, ...
    'arrivalInclination',     options.arrivalInclination, ...
    'arrivalArgOfPeriapsis',  options.arrivalArgOfPeriapsis, ...
    'tof',                    tof);
end

function [dv_loi, dv_plane, dv_trim] = biEllipticCapture( ...
        Moon, v_hyp_peri, r_peri, r_bi_apo, r_apo_final, inc_deg)
%BIELLIPTICCAPTURE Three-burn bi-elliptic lunar orbit capture.
%   Burn 1 (perilune):  hyperbola -> r_peri x r_bi_apo equatorial, no plane change
%   Burn 2 (apolune):   pure inc_deg plane change at r_bi_apo
%   Burn 3 (perilune):  adjust apoapsis to r_apo_final, no plane change
%
%   For r_bi_apo >= r_apo_final the v_bi_peri terms cancel and
%   dv_total = (v_hyp_peri - v_final_peri) + 2*v_bi_apo*sin(inc/2).

a_bi      = (r_peri + r_bi_apo) / 2;
v_bi_peri = sqrt(2*Moon.mu/r_peri - Moon.mu/a_bi);
v_bi_apo  = sqrt(2*Moon.mu/r_bi_apo - Moon.mu/a_bi);

a_final      = (r_peri + r_apo_final) / 2;
v_final_peri = sqrt(2*Moon.mu/r_peri - Moon.mu/a_final);

dv_loi   = v_hyp_peri - v_bi_peri;               % retrograde, no plane change
dv_plane = 2 * v_bi_apo * sin(deg2rad(inc_deg)/2); % pure plane change
dv_trim  = abs(v_bi_peri - v_final_peri);          % apoapsis adjustment
end

function dv = combinedBurnDV(v_before, v_after, inc_deg)
%COMBINEDBURNDV Delta-v for a simultaneous speed change and plane change.
%   Uses the vector triangle: dv = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(inc))
%   At inc = 0 this reduces to |v_after - v_before|.
dv = sqrt(v_before^2 + v_after^2 - 2*v_before*v_after*cos(deg2rad(inc_deg)));
end
