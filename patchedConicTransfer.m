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

% Quick path for Earth->Moon lunar transfer
if strcmpi(departBody.name, 'Earth') && strcmpi(arrivalBody.name, 'Moon')
    result = lunarTransfer(options.departureAltitude, options.arrivalAltitude, options);
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
% INTERPLANETARYTRANSFER simple Hohmann-style transfer between planets

% Assume circular coplanar orbits
r1 = departBody.a;
r2 = arrivalBody.a;

muSun = constants().Sun.mu;

% Transfer orbit
a_t = (r1 + r2)/2;

% Velocity in circular orbits
v1 = sqrt(muSun/r1);
v2 = sqrt(muSun/r2);

% Velocity at perihelion/aphelion of transfer orbit
v_peri = sqrt(2*muSun/r1 - muSun/a_t);
v_apo  = sqrt(2*muSun/r2 - muSun/a_t);

% Heliocentric delta-vs — these equal the hyperbolic excess speeds
% relative to each body for a coplanar Hohmann transfer.
v_inf_dep = abs(v_peri - v1);
v_inf_arr = abs(v2 - v_apo);

% Body-centric departure burn: combined TLI + plane change.
% departureInclination is the angle between the parking orbit plane
% and the ecliptic transfer plane.
r_park = departBody.radius + options.departureAltitude;
v_park = sqrt(departBody.mu / r_park);
v_hyp_dep = sqrt(v_inf_dep^2 + 2*departBody.mu/r_park);
dv_dep = combinedBurnDV(v_park, v_hyp_dep, options.departureInclination);

% Body-centric arrival burn: combined capture + plane change.
% arrivalInclination is the angle between the capture orbit plane
% and the arrival hyperbola plane.
r_cap = arrivalBody.radius + options.arrivalAltitude;
v_cap = sqrt(arrivalBody.mu / r_cap);
v_hyp_arr = sqrt(v_inf_arr^2 + 2*arrivalBody.mu/r_cap);
dv_arr = combinedBurnDV(v_hyp_arr, v_cap, options.arrivalInclination);

% Time of flight (half period)
tof = pi * sqrt(a_t^3/muSun);

% Phase angle (approx for circular coplanar)
phase_rad = pi*(1 - sqrt((2*r2)/(r1 + r2)));
phase_deg = rad2deg(phase_rad);

res = struct();
res.deltaV = dv_dep + dv_arr;
res.tof = tof;
res.phaseAngle = phase_deg;
res.details = struct('dvDeparture', dv_dep, ...
                         'dvArrival',   dv_arr, ...
                         'tof', tof, ...
                         'r1', r1, ...
                         'r2', r2, ...
                         'aTransfer',   a_t, ...
                         'vInfDepart',  v_inf_dep, ...
                         'vInfArrive',  v_inf_arr, ...
                         'rParkDepart', r_park, ...
                         'rParkArrive', r_cap, ...
                         'departureInclination', options.departureInclination, ...
                         'arrivalInclination',   options.arrivalInclination);
end

function dv = combinedBurnDV(v_before, v_after, inc_deg)
%COMBINEDBURNDV Delta-v for a simultaneous speed change and plane change.
%   Uses the vector triangle: dv = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(inc))
%   At inc = 0 this reduces to |v_after - v_before|.
dv = sqrt(v_before^2 + v_after^2 - 2*v_before*v_after*cos(deg2rad(inc_deg)));
end
