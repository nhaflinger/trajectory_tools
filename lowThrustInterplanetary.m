function result = lowThrustInterplanetary(departBody, arrivalBody, departJD, tofDays, options)
%LOWTHRUSTINTERPLANETARY  Three-phase low-thrust interplanetary transfer.
%
%   result = lowThrustInterplanetary(departBody, arrivalBody, departJD, tofDays)
%   result = lowThrustInterplanetary(departBody, arrivalBody, departJD, tofDays, options)
%
%   Models a complete low-thrust mission as three sequential Edelbaum spirals:
%
%     Phase 1 — Departure escape spiral
%       Parking orbit (r_park_dep) up to departure body's sphere of influence.
%       Central body: departBody.
%
%     Phase 2 — Heliocentric cruise
%       From departure body's heliocentric distance to arrival body's distance.
%       Uses actual planet positions at departure/arrival dates (not just
%       semi-major axes), so eccentricity and phasing affect the result.
%       Central body: Sun.
%
%     Phase 3 — Arrival capture spiral
%       From arrival body's sphere of influence down to parking orbit.
%       Central body: arrivalBody.
%
%   Mass is propagated sequentially: each phase starts with the dry+residual
%   mass left after the previous phase.  The Edelbaum formula assumes the
%   transfer takes many orbital revolutions (optimal slow spiral); it is
%   accurate for low-thrust, high-Isp systems and gives an optimistic lower
%   bound for higher-thrust cases.
%
%   NOTE: Edelbaum for the heliocentric phase yields |v1 - v2|, independent
%   of TOF at fixed planet positions.  TOF sensitivity comes mainly from
%   the changing planet distances with the departure and arrival dates.
%   Use porkChopLowThrust() to survey the full (departureDate, TOF) space.
%
%   Inputs:
%     departBody  - body struct from constants() — departure planet
%     arrivalBody - body struct from constants() — destination planet
%     departJD    - departure Julian Date
%     tofDays     - total mission time of flight (days)
%     options     (optional struct):
%       .thrustN               engine thrust (N)                      [0]
%       .isp                   specific impulse (s)                  [3000]
%       .wetMass               initial spacecraft wet mass (kg)      [1000]
%       .departureAltitude     parking orbit altitude at depart (km)  [200]
%       .arrivalAltitude       parking orbit altitude at arrival (km) [400]
%       .nStepsPerOrbit        RK4 steps per orbit (per phase)         [50]
%       .nOutputPoints         trajectory output points (per phase)  [1000]
%       .skipDepartureSpiral   true => launch vehicle provides escape  [false]
%                              Phase 1 ΔV set to 0; mass unchanged
%       .skipArrivalSpiral     true => no capture into parking orbit   [false]
%                              Phase 3 ΔV set to 0; mass unchanged
%
%   Outputs (result struct):
%     .deltaV          total mission ΔV: sum of three phases (km/s)
%     .tof             total TOF (s)      [NaN if thrustN = 0]
%     .tofDays         total TOF (days)
%     .propellantMass  total propellant consumed (kg) [NaN if thrustN = 0]
%     .finalMass       final spacecraft mass (kg)
%     .details         per-phase breakdown (dvDeparture/Heliocentric/Arrival,
%                      tofDays*, propellant*)
%     .phases          struct with full lowThrustSpiral results for each phase
%                      (.departure  .heliocentric  .arrival)

if nargin < 5, options = struct(); end
if ~isfield(options,'thrustN'),            options.thrustN            = 0;    end
if ~isfield(options,'isp'),                options.isp                = 3000; end
if ~isfield(options,'wetMass'),            options.wetMass            = 1000; end
if ~isfield(options,'departureAltitude'),  options.departureAltitude  = 200;  end
if ~isfield(options,'arrivalAltitude'),    options.arrivalAltitude    = 400;  end
if ~isfield(options,'nStepsPerOrbit'),     options.nStepsPerOrbit     = 50;   end
if ~isfield(options,'nOutputPoints'),      options.nOutputPoints      = 1000; end
if ~isfield(options,'skipDepartureSpiral'), options.skipDepartureSpiral = false; end
if ~isfield(options,'skipArrivalSpiral'),   options.skipArrivalSpiral   = false; end

consts = constants();
muSun  = consts.Sun.mu;

% ---- Sphere-of-influence radii (same formula as constants.m) --------
r_SOI_dep = departBody.a * (departBody.mu / muSun)^(2/5);
r_SOI_arr = arrivalBody.a * (arrivalBody.mu / muSun)^(2/5);

% ---- Parking orbit radii --------------------------------------------
r_park_dep = departBody.radius + options.departureAltitude;
r_park_arr = arrivalBody.radius + options.arrivalAltitude;

% ---- Actual heliocentric distances at departure / arrival -----------
% Using true eccentric orbits so planet eccentricity (e.g. Mars e=0.093)
% contributes to the heliocentric speed difference.
arrJD = departJD + tofDays;
try
    r1_vec = orbitalState(departBody, departJD);
    r2_vec = orbitalState(arrivalBody, arrJD);
    r1_helio = norm(r1_vec);
    r2_helio = norm(r2_vec);
catch
    % Fall back to semi-major axis if ephemeris fails
    r1_helio = departBody.a;
    r2_helio = arrivalBody.a;
end

% ---- Build thruster options for each phase (mass threaded through) --
g0 = 9.80665e-3;   % km/s^2

sunBody = consts.Sun;

function o = phaseOpts(wetMass_kg)
    o               = struct();
    o.thrustN       = options.thrustN;
    o.isp           = options.isp;
    o.wetMass       = wetMass_kg;
    o.nStepsPerOrbit = options.nStepsPerOrbit;
    o.nOutputPoints  = options.nOutputPoints;
end

% ---- Phase 1: departure body spiral (parking -> SOI) ----------------
if options.skipDepartureSpiral
    phase1 = struct('deltaV', 0, 'tof', 0, 'tofDays', 0, ...
                    'propellantMass', 0, 'finalMass', options.wetMass, ...
                    'trajectory', struct());
    m_after1 = options.wetMass;
else
    opts1   = phaseOpts(options.wetMass);
    phase1  = lowThrustSpiral(departBody, r_park_dep, r_SOI_dep, opts1);

    m_after1 = options.wetMass;
    if options.thrustN > 0 && ~isnan(phase1.finalMass)
        m_after1 = phase1.finalMass;
    else
        F_kN = options.thrustN * 1e-3;
        if F_kN > 0
            m_after1 = options.wetMass * exp(-phase1.deltaV / (options.isp * g0));
        end
    end
end

% ---- Phase 2: heliocentric spiral (dep orbit -> arr orbit) ----------
opts2  = phaseOpts(m_after1);
phase2 = lowThrustSpiral(sunBody, r1_helio, r2_helio, opts2);

m_after2 = m_after1;
if options.thrustN > 0 && ~isnan(phase2.finalMass)
    m_after2 = phase2.finalMass;
else
    F_kN = options.thrustN * 1e-3;
    if F_kN > 0
        m_after2 = m_after1 * exp(-phase2.deltaV / (options.isp * g0));
    end
end

% ---- Phase 3: arrival body spiral (SOI -> parking) ------------------
if options.skipArrivalSpiral
    phase3 = struct('deltaV', 0, 'tof', 0, 'tofDays', 0, ...
                    'propellantMass', 0, 'finalMass', m_after2, ...
                    'trajectory', struct());
    m_final = m_after2;
else
    opts3  = phaseOpts(m_after2);
    phase3 = lowThrustSpiral(arrivalBody, r_SOI_arr, r_park_arr, opts3);

    m_final = m_after2;
    if options.thrustN > 0 && ~isnan(phase3.finalMass)
        m_final = phase3.finalMass;
    else
        F_kN = options.thrustN * 1e-3;
        if F_kN > 0
            m_final = m_after2 * exp(-phase3.deltaV / (options.isp * g0));
        end
    end
end

% ---- Totals ---------------------------------------------------------
dv_total  = phase1.deltaV + phase2.deltaV + phase3.deltaV;
prop_total = options.wetMass - m_final;

if options.thrustN > 0
    tofs = [phase1.tof, phase2.tof, phase3.tof];
    tof_total = sum(tofs(isfinite(tofs)));
else
    tof_total = NaN;
end

% ---- Pack result ----------------------------------------------------
result             = struct();
result.deltaV      = dv_total;
result.tof         = tof_total;
result.tofDays     = tof_total / 86400;
result.propellantMass = prop_total;
result.finalMass   = m_final;

result.details = struct( ...
    'dvDeparture',        phase1.deltaV,   ...
    'dvHeliocentric',     phase2.deltaV,   ...
    'dvArrival',          phase3.deltaV,   ...
    'tofDaysDeparture',   phase1.tofDays,  ...
    'tofDaysHeliocentric', phase2.tofDays, ...
    'tofDaysArrival',     phase3.tofDays,  ...
    'propDeparture',      phase1.propellantMass,    ...
    'propHeliocentric',   phase2.propellantMass,    ...
    'propArrival',        phase3.propellantMass,    ...
    'r1_helio',           r1_helio,        ...
    'r2_helio',           r2_helio,        ...
    'r_SOI_dep',          r_SOI_dep,       ...
    'r_SOI_arr',          r_SOI_arr,       ...
    'r_park_dep',         r_park_dep,      ...
    'r_park_arr',         r_park_arr,      ...
    'departJD',           departJD,        ...
    'arrJD',              arrJD);

result.phases = struct( ...
    'departure',    phase1,  ...
    'heliocentric', phase2,  ...
    'arrival',      phase3);
end
