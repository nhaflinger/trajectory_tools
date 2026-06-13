# API Reference

This is a concise, function-level API reference for the trajectory_tools MATLAB library — call signatures, inputs, and outputs for each core function. For tutorials, worked examples, and theory/background, see [README.md](../README.md).

## Constants & Utilities

### `SSBodies = constants()`
Returns a struct of physical/orbital constants for the Sun, planets, and Moon, for use in patched-conic and interplanetary calculations. Elements (a, e, i, Omega, omega_peri, M0) are J2000 Keplerian elements (Standish 1992, valid 1800–2050).

**Outputs**
- `SSBodies.Constants` — `G`, `AU` (km), `day` (s)
- `SSBodies.Sun`, `.Earth`, `.Moon`, `.Mercury`, ... — each with `name`, `mu` (km^3/s^2), `radius` (km), and (for orbiting bodies) `a`, `e`, `inclination`, `Omega`, `omega_peri`, `M0`, `obliquity`

### `jd = julianDate(year, month, day, hour, minute, second)`
Converts a calendar date/time (UTC) to Julian Date using the standard Gregorian-calendar algorithm. `hour`, `minute`, `second` default to 0.

**Inputs**
- `year, month, day` — calendar date (month 1-12)
- `hour, minute, second` — optional time of day

**Outputs**
- `jd` — Julian Date (scalar)

### `[r_sun_hat, r_sun_km] = sunPosition(jd)`
Computes the Earth-to-Sun direction in ECI J2000 using a low-precision solar almanac (~0.01° accuracy; Vallado Algorithm 29). Sufficient for beta-angle and eclipse calculations.

**Inputs**
- `jd` — Julian Date (scalar)

**Outputs**
- `r_sun_hat` — unit vector Earth→Sun, ECI J2000, 3x1
- `r_sun_km` — full Earth→Sun position vector (km), 3x1

### `varargout = wgs84Geodetic(varargin)`
Converts between ECEF Cartesian coordinates and WGS84 geodetic latitude/longitude/altitude, using the Bowring iterative method (ECEF→geodetic) or a closed-form expression (geodetic→ECEF).

**Usage**
- `[lat_deg, lon_deg, alt_km] = wgs84Geodetic(r_ecef_km)` — `r_ecef_km` is 3x1/1x3, or 3xN for vectorized output
- `r_ecef_km = wgs84Geodetic(lat_deg, lon_deg, alt_km, 'inverse')`

---

## Orbital State, Elements & Propagation

### `[r_vec, v_vec] = orbitalState(body, jd, epoch0)`
Computes a heliocentric ecliptic J2000 state vector for a solar-system body by propagating its J2000 Keplerian elements to `jd`. Handles elliptical, circular-coplanar fallback (if `Omega`/`omega_peri`/`M0` absent), and hyperbolic (`e >= 1`, requires `body.t_peri_jd`) orbits.

**Inputs**
- `body` — struct from `constants()` with `a, e, inclination, Omega, omega_peri, M0` (or `t_peri_jd` for hyperbolic)
- `jd` — Julian Date to evaluate
- `epoch0` — optional reference epoch JD (default J2000 = 2451545.0)

**Outputs**
- `r_vec, v_vec` — heliocentric ecliptic position (km) and velocity (km/s), 3x1

### `[r, v] = orbitalStateCircular(body, epoch, epoch0)`
Thin wrapper/delegate to `orbitalState` — kept for backward compatibility. Produces an eccentric 3-D solution when full elements are present, else falls back to circular-coplanar.

**Inputs**
- `body, epoch, epoch0` — same as `orbitalState`

**Outputs**
- `r, v` — heliocentric position (km) and velocity (km/s), 3x1

### `[r_vec, v_vec] = coe2eci(a, e, i, RAAN, omega, nu)`
Converts classical (Keplerian) orbital elements to an Earth-centered inertial (ECI) state vector via the perifocal (PQW) frame and standard 3-1-3 rotation. Uses Earth `mu = 398600.4418 km^3/s^2`.

**Inputs**
- `a` — semi-major axis (km), `e` — eccentricity
- `i, RAAN, omega, nu` — inclination, RAAN, argument of perigee, true anomaly (all deg)

**Outputs**
- `r_vec` — ECI position (km, 3x1), `v_vec` — ECI velocity (km/s, 3x1)

### `coe = eci2coe(r_vec, v_vec)`
Converts an ECI state vector to classical orbital elements for an Earth orbit.

**Inputs**
- `r_vec, v_vec` — ECI position (km) and velocity (km/s), 3x1 or 1x3

**Outputs (struct fields)**
- `a` (km), `e`, `i` (deg), `RAAN` (deg), `omega` (deg), `nu` (deg), `M` (deg), `p` (semi-latus rectum, km), `h_vec` (angular momentum vector, km^2/s)

### `orb = earthOrbit(type, varargin)`
Factory function that builds a standardized Earth-orbit definition struct from several input conventions (COE, ECI state, circular, sun-synchronous, GEO, Molniya, or TLE).

**Usage forms**
- `earthOrbit('coe', a_km, e, i_deg, RAAN_deg, omega_deg, M0_deg [, epoch_jd])`
- `earthOrbit('eci', r_vec, v_vec [, epoch_jd])`
- `earthOrbit('circular', alt_km, i_deg [, RAAN_deg])`
- `earthOrbit('sso', alt_km [, e, omega_deg, ltan_hrs, epoch_jd])`
- `earthOrbit('geo' [, lon_deg])`, `earthOrbit('molniya', RAAN_deg [, omega_deg])`, `earthOrbit('tle', line1, line2)`

**Key outputs (struct fields)**
- `a, e, i, RAAN, omega, M0, epoch_jd, period` (s)
- `alt_peri, alt_apo` (km), `r_vec, v_vec` (ECI at epoch)
- `type`, `ltan_hrs`, `sun_sync_RAAN` (NaN if not applicable)

### `result = propagateOrbit(orb, duration_s, varargin)`
Propagates an `earthOrbit` struct forward in time and returns a time-history struct, with selectable fidelity (analytical Kepler/J2 secular, or numerical with J2/J3/J4/drag/SRP).

**Inputs**
- `orb` — orbit struct from `earthOrbit()`
- `duration_s` — propagation duration (s)
- Options: `'Method'` — `'kepler'|'j2'(default)|'numerical'|'drag'|'j3'|'j4'|'srp'|'full'`; `'StepSize'` (s); `'CdAm'`, `'CR_Am'`; `'OutputCOE'` (bool)

**Outputs (struct fields)**
- `t` (s, Nx1), `r_eci`, `v_eci` (3xN), `lat, lon, alt` (deg/deg/km, Nx1), `orb`
- If `OutputCOE=true`: `a, e, i_deg, RAAN, omega, M` (each Nx1)

---

## Lambert Solver & Patched-Conic Transfers

### `[v1, v2] = lambertSolver(r1_vec, r2_vec, tof, mu, isLongWay)`
Universal-variable Lambert solver (Bate-Mueller-White §5.3) using a bisection search on the universal variable `z`, since TOF(z) is monotonic. Returns the departure/arrival velocities for the transfer orbit connecting two position vectors in a given time of flight.

**Inputs**
- `r1_vec, r2_vec` — departure/arrival position vectors (km, 3x1)
- `tof` — time of flight (s)
- `mu` — gravitational parameter of central body (km^3/s^2)
- `isLongWay` — true for >180° transfer angle (default false)

**Outputs**
- `v1, v2` — departure/arrival velocity vectors (km/s, 3x1)

### `result = patchedConicTransfer(departBody, arrivalBody, options)`
Simplified patched-conic transfer between two solar-system bodies, assuming circular coplanar orbits with a Hohmann-like Lambert leg. For Earth→Moon, uses a two-phase Earth-centered model (parking orbit → TLI → lunar capture).

**Inputs (options)**
- `departureAltitude`, `arrivalAltitude` (km, default 200)
- `departureInclination`, `arrivalInclination` (deg, default 0)
- `arrivalApogeeAltitude` (km, default = arrivalAltitude — elliptical capture orbit)
- `nRevs` — Lambert solver revolutions (default 0)

**Outputs (struct fields)**
- `deltaV` (km/s, total), `tof` (s), `phaseAngle` (deg), `details` (breakdown struct)

### `plotPatchedConic(result, departBody, arrivalBody)`
Visualizes a `patchedConicTransfer` result: a bar chart of ΔV components plus a 2-D (or 3-D for Earth-Moon) view of the transfer trajectory.

**Inputs**
- `result` — output of `patchedConicTransfer`
- `departBody, arrivalBody` — body structs from `constants()`

### `porkChopPlot(departBody, arrivalBody, departJD, tofDays, options)`
Generates a pork-chop contour plot (departure date vs. time-of-flight, colored by ΔV) using Lambert solutions on circular orbits for each grid point.

**Inputs**
- `departJD` — vector of departure Julian Dates
- `tofDays` — vector of time-of-flight values (days)
- `options.epoch0` (default J2000), `options.mu` (default Sun mu), `options.fixPlane` (default true)

### `best = findBestLaunchDate(departBody, arrivalBody, jdStart, jdEnd, tofDays, options)`
Grid-searches departure dates and TOFs (Lambert on circular orbits) to find the minimum-ΔV transfer in a window. If `jdStart`/`jdEnd` are scalars, uses 5-day steps.

**Inputs**
- `jdStart, jdEnd` — departure-date window (JD), or a full `departJD` vector
- `tofDays` — vector of TOF values (days)
- `options.epoch0`, `options.mu` — as above

**Outputs (struct fields)**
- `departureJD`, `tofDays`, `deltaV` (minimum), `dvGrid` (full matrix)

---

## Gravity Assist & Multi-Body Trajectories

### `result = gravityAssist(body, v_inf_in_vec, v_inf_out_vec, options)`
Analyzes a single gravity-assist flyby from incoming/outgoing heliocentric v∞ vectors. Determines deflection angle, required periapsis radius, feasibility, and — if `|v∞_in| ≠ |v∞_out|` — the powered-flyby ΔV (Oberth burn at periapsis).

**Inputs**
- `body` — struct from `constants()` (`.mu`, `.radius`, `.name`)
- `v_inf_in_vec, v_inf_out_vec` — heliocentric v∞ vectors (km/s, 3x1)
- `options.atmosphereAltitude` (km, default 0), `options.flybyAltitude` (km, default 300)

**Outputs (struct fields)**
- `v_inf_in, v_inf_out` (km/s), `deflection` (deg), `r_periapsis` (km), `altitude` (km)
- `isFeasible`, `isPowered`, `dvPowered` (km/s), `maxDeflection` (deg)
- `v_inf_in_vec, v_inf_out_vec`

### `result = flybySequence(bodies, departureJD, tofDays, options)`
Chains N−1 Lambert legs through N bodies (departure, intermediate flybys, arrival), analyzing flyby geometry/feasibility at each intermediate body via `gravityAssist`.

**Inputs**
- `bodies` — 1xN cell array of body structs from `constants()`
- `departureJD` — departure JD from the first body
- `tofDays` — 1x(N-1) per-leg time of flight (days)
- `options`: `departureAltitude` (200), `arrivalAltitude` (400), `arrivalApogeeAltitude`, `departureInclination`/`arrivalInclination` (0), `flybyAltitudes` (300 each), `atmosphereAltitudes` (0 each), `transferTypes` (cell of `'type1'|'type2'`, default `'type1'`)

**Outputs (struct fields)**
- `legs` — 1x(N-1) per-leg struct array (`r1_vec, v1_body, v1_transfer, v_inf_depart`, etc.)
- `flybys` — 1x(N-2) per-flyby struct array (see `gravityAssist`)
- `deltaV` (total, km/s), `deltaVBurns` (no TCM), `tof` (days), `details` (budget breakdown)

### `best = findBestFlybyWindow(bodies, jdStart, jdEnd, tofRanges, options)`
Grid-searches departure dates and per-leg TOFs to minimize a heliocentric ΔV proxy (v∞ departure + Σ powered-flyby ΔV + v∞ arrival).

**Inputs**
- `bodies` — 1xN cell array (departure, flybys, arrival)
- `jdStart, jdEnd` — departure date search window
- `tofRanges` — (N-1)x2 `[min, max]` days per leg, or explicit (N-1)xM TOF grid
- `options.nDepDates` (50), `options.nTofPoints` (25), `options.flybyAltitudes`, `options.transferTypes`

**Outputs (struct fields)**
- `departureJD`, `tofDays` (1x(N-1)), `deltaVProxy` (km/s), `departJD`, `tofGrids`

### `plotFlybySequence(result, bodies)`
3-D heliocentric visualization of a gravity-assist sequence: ecliptic plane, Sun, full body orbits, active Lambert arcs per leg, encounter markers, and a ΔV/flyby summary annotation.

**Inputs**
- `result` — output of `flybySequence()`
- `bodies` — 1xN cell array matching the order passed to `flybySequence`

### `[fig, ax] = tisserandGraph(bodies, options)`
Plots a Tisserand graph (iso-v∞ contours per body on the periapsis-vs-apoapsis plane) for gravity-assist mission design. Free gravity assists move the spacecraft along a body's contour; deep-space maneuvers move it between contours.

**Inputs**
- `bodies` — 1xN cell array of body structs
- `options.vInfValues` (km/s, default `[0.5 1 2 3 5 10 20]`), `options.nGrid` (400), `options.rpLim`/`raLim` (AU, auto), `options.trajectory` (overlay a `flybySequence` result), `options.showLabels` (true)

**Outputs**
- `fig` — figure handle (and axes `ax`)

### `result = resonantOrbits(body, varargin)`
Computes p:q resonant-orbit parameters (spacecraft completes p revs while the planet completes q) for gravity-assist sequencing, including required semi-major axis, period, minimum tangent-encounter v∞, and apsis distances.

**Inputs**
- `body` — flyby planet struct from `constants()`
- Options: `maxN` (max p/q integer, default 5), `ax` (axes handle to annotate a Tisserand graph), `print` (default true)

**Outputs**
- `result` — struct array (one per resonance) with fields `label` (`'p:q'`), `p, q`, `a_sc` (km), `T_days`, `v_inf` (km/s), `r_p, r_a` (km)

### `fig = porkChopSequence(bodies, jdStart, jdEnd, tofRanges, options)`
Generates one pork-chop subplot per leg of a multi-body flyby sequence, showing arrival v∞ vs. departure date and TOF, to identify departure windows where consecutive legs' v∞ values match (free gravity assist).

**Inputs**
- `bodies` — 1xN cell array
- `jdStart, jdEnd` — leg-1 departure window
- `tofRanges` — (N-1)x2 `[min, max]` days per leg
- `options.nDep`/`nTof` (60), `options.transferTypes`, `options.markBest` (struct with `.departureJD`/`.tofDays`), `options.cLimPct` (`[5 85]`)

**Outputs**
- `fig` — figure handle

---

## Low-Thrust Transfers

### `result = lowThrustSpiral(body, r0_km, r1_km, options)`
Computes ΔV and propellant budget for a continuous low-thrust transfer between two circular orbits, combining the Edelbaum (1961) analytical ΔV with an optional tangential-thrust RK4 propagation and Tsiolkovsky mass accounting.

**Inputs**
- `body` — central body struct from `constants()`
- `r0_km, r1_km` — initial/target circular orbit radii (km from body center)
- `options`: `inclinationChangeDeg` (0), `thrustN` (0 → Edelbaum-only), `isp` (3000 s), `wetMass` (1000 kg), `nStepsPerOrbit` (50), `nOutputPoints` (2000)

**Outputs (struct fields)**
- `deltaV` (km/s, Edelbaum), `tof`/`tofDays` (NaN if `thrustN=0`)
- `propellantMass`, `finalMass` (kg, NaN if `thrustN=0`)
- `details`, `trajectory` (`.t .x .y .vx .vy .r .speed .mass`, empty if `thrustN=0`)

### `fig = plotLowThrustSpiral(result, body, options)`
Three-panel visualization of a low-thrust spiral: 2-D spiral path (colored by time), altitude-vs-time, and mass-vs-time.

**Inputs**
- `result` — output of `lowThrustSpiral`
- `body` — central body struct
- `options.title`, `options.units` (`'km'` default or `'AU'`)

**Outputs**
- `fig` — figure handle (empty if no trajectory was computed)

### `result = lowThrustInterplanetary(departBody, arrivalBody, departJD, tofDays, options)`
Models a full low-thrust interplanetary mission as three sequential Edelbaum spirals: departure escape spiral, heliocentric cruise (using actual planet positions), and arrival capture spiral. Mass is propagated sequentially across phases.

**Inputs**
- `departBody, arrivalBody` — body structs from `constants()`
- `departJD` — departure Julian Date, `tofDays` — total mission TOF (days)
- `options`: `thrustN` (0), `isp` (3000), `wetMass` (1000 kg), `departureAltitude`/`arrivalAltitude` (200/400 km), `nStepsPerOrbit` (50), `nOutputPoints` (1000), `skipDepartureSpiral`/`skipArrivalSpiral` (false — assume LV escape / no capture)

**Outputs (struct fields)**
- `deltaV` (total, km/s), `tof`/`tofDays`, `propellantMass`, `finalMass` (kg)
- `details` — per-phase ΔV/TOF/propellant breakdown
- `phases` — full `lowThrustSpiral` results per phase

### `fig = plotLowThrustInterplanetary(result, departBody, arrivalBody, options)`
Four-panel visualization: heliocentric cruise view (planet orbits + colored trajectory), departure escape spiral, arrival capture spiral, and mass-vs-time across all phases.

**Inputs**
- `result` — output of `lowThrustInterplanetary`
- `departBody, arrivalBody` — body structs
- `options.AU`, `options.orbitTraceN` (360), `options.lvEscapeDV` (NaN), `options.lvParkAlt` (200 km)

**Outputs**
- `fig` — figure handle

### `fig = porkChopLowThrust(departBody, arrivalBody, departJD, tofDays, options)`
Low-thrust pork-chop plot: for each (departure date, TOF) grid point, computes the total three-phase Edelbaum mission ΔV/propellant fraction using actual planet positions, with an optional impulsive-Lambert contour overlay for comparison.

**Inputs**
- `departJD` — vector of departure JDs, `tofDays` — vector of TOF values (days)
- `options`: `isp` (3000), `wetMass` (1000 kg), `departureAltitude`/`arrivalAltitude` (200/400), `colorMode` (`'propFrac'` default or `'deltaV'`), `cLimPct` (`[5 85]`), `showImpulsive` (true), `impulsiveContours`

**Outputs**
- `fig` — figure handle

---

## Ground Tracks, Coverage & Access

### `fig = plotGroundTrack(orb, varargin)`
Mercator ground-track plot for an Earth orbit, with optional J2 RAAN/argument-of-perigee drift, multi-orbit coloring, ground station markers, and ascending-node annotations.

**Inputs**
- `orb` — struct from `earthOrbit()`
- Options: `'NumOrbits'` (3), `'StepsPerOrbit'` (360), `'J2'` (true), `'ColorByOrbit'` (true), `'GroundStations'` (struct array with `.lat .lon .name`), `'ShowNodes'` (true)

**Outputs**
- `fig` — figure handle

### `fig = plotOrbit3D(orb, varargin)`
3-D ECI visualization of an Earth orbit (equatorial plane, node line, perigee marker, optional J2 drift over multiple revolutions).

**Inputs**
- `orb` — struct from `earthOrbit()`
- Options: `'NumOrbits'` (1), `'StepsPerOrbit'` (720), `'J2'` (false), `'ShowEquator'` (true), `'ShowNodeLine'` (true), `'ShowPerigee'` (true)

**Outputs**
- `fig` — figure handle

### `downloadCoastlines()`
One-time helper that downloads the Natural Earth 110m coastline GeoJSON and caches it as `coastlines_cache.mat` (used by `plotGroundTrack`). Run once before plotting ground tracks, or it is invoked automatically.

### `cov = coverageAnalysis(orb, duration_s, varargin)`
Propagates an orbit and computes coverage fraction and revisit statistics over a lat/lon grid, based on a minimum elevation-angle visibility criterion.

**Inputs**
- `orb` — struct from `earthOrbit()`, `duration_s` — analysis duration (s)
- Options: `'GridRes'` (5°), `'LatLim'`/`'LonLim'` (`[-90 90]`/`[-180 180]`), `'MinElev'` (5°), `'Method'` (`'j2'`), `'StepSize'` (auto)

**Outputs (struct fields)**
- `lat_vec, lon_vec` — grid vectors (deg)
- `coverage_frac`, `revisit_mean_hr`, `revisit_max_hr`, `n_passes` — MxN matrices
- `duration_hr`, `orb`

### `fig = plotCoverage(cov, varargin)`
Plots `coverageAnalysis` (or `constellationCoverage`) results on a world map, either coverage fraction or revisit time.

**Inputs**
- `cov` — struct from `coverageAnalysis` or `constellationCoverage`
- Options: `'Quantity'` (`'coverage'` default or `'revisit'`), `'GroundStations'` (struct array)

**Outputs**
- `fig` — figure handle

### `wins = accessWindows(orb, target_lat, target_lon, duration_s, varargin)`
Propagates an orbit and computes discrete visibility windows over a fixed ground target based on a minimum elevation angle.

**Inputs**
- `orb`, `target_lat, target_lon` (deg), `duration_s` (s)
- Options: `'MinElev'` (5°), `'Method'` (`'j2'`), `'StepSize'` (auto)

**Outputs**
- `wins` — struct array with `start_s, stop_s, duration_s, max_elev_deg, start_jd, stop_jd`

### `[az, el, range_km] = topocentricAzEl(obs_lat, obs_lon, obs_alt_km, r_sat_eci, jd)`
Computes azimuth, elevation, and slant range from a ground observer to a satellite by rotating the satellite ECI position to ECEF (via GMST) and projecting into the observer's ENU frame.

**Inputs**
- `obs_lat, obs_lon` (deg), `obs_alt_km` (km)
- `r_sat_eci` — satellite ECI position (km), 3x1 or 3xN
- `jd` — Julian Date, scalar or 1xN

**Outputs**
- `az, el` (deg, 1xN), `range_km` (1xN)

### `[fig, summary] = plotAccessWindowGantt(orb, ground_stations, duration_hr, varargin)`
Gantt-style timeline of access windows across multiple ground stations for a given orbit.

**Inputs**
- `orb` — struct from `earthOrbit()`
- `ground_stations` — struct array with `.lat, .lon, .alt (optional), .name`
- `duration_hr` — analysis duration (hours)
- Options: `'MinElevation'` (10°), `'StepSize'` (30 s), `'ColorByElevation'` (true), `'ShowMaxEl'` (true), `'TimeFormat'` (`'hours'`|`'minutes'`)

**Outputs**
- `fig` — figure handle
- `summary` — struct array per station: `name, n_passes, total_contact_min, mean_duration_min, max_elevation_deg, windows`

---

## Orbital Maneuvers

### `res = hohmannTransfer(arg1, arg2, varargin)`
Computes a Hohmann (or bi-elliptic) transfer between two circular Earth orbits. Each argument may be an altitude (km) or an `earthOrbit` struct.

**Inputs**
- `arg1, arg2` — altitude (km) or `earthOrbit` struct
- `'Bielliptic', r_int_km` — optional, intermediate apoapsis radius (km) for a bi-elliptic transfer

**Outputs (struct fields)**
- `dv1_km_s, dv2_km_s, dv_total_km_s`, `tof_s, tof_min`
- `alt1_km, alt2_km, r1_km, r2_km`, `a_transfer_km` (or `a1_km, a2_km` for bi-elliptic)
- `type` (`'hohmann'`|`'bielliptic'`), plus `dv3_km_s` for bi-elliptic

### `res = planeChange(orb, delta_i_deg, varargin)`
Computes the ΔV for a pure inclination-change maneuver, or a combined plane-change + Hohmann altitude change at the transfer apogee.

**Inputs**
- `orb` — `earthOrbit` struct (circular orbit assumed)
- `delta_i_deg` — inclination change magnitude (deg, >= 0)
- `'Type'` — `'pure'` (default) or `'combined'`; `'AltFinal_km'` required for `'combined'`

**Outputs (struct fields)**
- `dv_km_s`, `dv_separate_km_s`/`savings_km_s` (combined only), `delta_i_deg`, `type`, `alt_km`

### `res = phasingManeuver(orb, delta_phase_deg, n_phasing_orbits)`
Computes the two-burn ΔV and timing for a phasing maneuver that shifts a spacecraft's position by `delta_phase_deg` over `n_phasing_orbits`.

**Inputs**
- `orb` — `earthOrbit` struct (circular)
- `delta_phase_deg` — desired phase change (deg); positive = catch up, negative = open gap
- `n_phasing_orbits` — number of phasing-orbit revolutions

**Outputs (struct fields)**
- `dv_each_km_s, dv_total_km_s`, `a_phase_km`, `alt_peri_phase_km`/`alt_apo_phase_km`
- `period_phase_s`/`period_phase_min`, `time_to_complete_s`/`_hr`, `direction` (`'catch_up'`|`'open_gap'`)

### `res = deorbitBurn(orb_or_alt, varargin)`
Computes the retrograde ΔV to lower perigee from a circular orbit to a target reentry altitude.

**Inputs**
- `orb_or_alt` — `earthOrbit` struct or scalar circular altitude (km)
- `'TargetPeriAlt_km'` — target perigee altitude after burn (default 80 km)

**Outputs (struct fields)**
- `dv_km_s`/`dv_m_s`, `alt_km`, `target_peri_alt_km`, `tof_to_peri_s`/`_min`, `a_deorbit_km`, `compliant_25yr` (true if alt <= 600 km)

### `wins = launchWindow(launch_lat, launch_lon, target_inc, target_RAAN, start_jd, varargin)`
Computes launch windows (times, azimuths, achieved RAAN/inclination) for reaching a target orbital plane from a ground site, accounting for optional target-RAAN drift.

**Inputs**
- `launch_lat, launch_lon` (deg), `target_inc, target_RAAN` (deg), `start_jd`
- Options: `'NDays'` (1), `'SearchStep_s'` (60), `'RAANDrift_degPerDay'` (0)

**Outputs**
- `wins` — struct array with `jd, date_str, type` (`'ascending'|'descending'`), `azimuth_deg, achieved_RAAN_deg, achieved_inc_deg, LMST_deg`

### `fig = plotLaunchWindow(wins, launch_lat, launch_lon, target_inc, target_RAAN, start_jd, n_days)`
Two-panel visualization of `launchWindow` results: LMST vs. time with required-LMST lines and window markers (top), launch azimuth per window (bottom).

**Inputs**
- `wins` — output of `launchWindow`
- Remaining args mirror `launchWindow`'s inputs plus `n_days` (search window length)

**Outputs**
- `fig` — figure handle

### `[result, fig] = dogLegTrade(launch_lat, launch_lon, target_inc, target_RAAN, start_jd, varargin)`
Sweeps launch times over a day and computes the dog-leg ΔV penalty (plane-change cost) for launching at non-optimal times relative to the zero-cost window: `dv_dogleg = 2*v_insertion*sin(delta_RAAN/2)`.

**Inputs**
- `launch_lat, launch_lon, target_inc, target_RAAN, start_jd` — as in `launchWindow`
- Options: `'NDays'` (1), `'InsertionAlt_km'` (300), `'Plot'` (false)

**Outputs (struct fields)**
- `t_hr, jd, dv_dogleg_m_s, RAAN_asc, RAAN_desc, window_times_hr, min_dv_m_s, max_dv_m_s`
- `fig` — figure handle (if `'Plot'` true)

### `budget = geoStationKeeping(varargin)`
Computes the annual North-South (inclination, due to lunar/solar perturbations) and East-West (longitude drift, due to triaxiality) station-keeping ΔV budget for a GEO satellite.

**Inputs**
- `lon_deg` — GEO longitude slot (deg E, default 0)
- `year` — reference year (default current year)

**Outputs (struct fields)**
- `dv_ns_km_s`/`_m_s`, `dv_ew_km_s`/`_m_s`, `dv_total_km_s`/`_m_s`
- `lon_deg`, `v_geo_km_s`, `a_geo_km`, `lifetime_yr_per_10ms`, `deadband_ew_deg`

---

## Relative Motion (Clohessy-Wiltshire)

### `[r_hist, v_hist] = cwPropagate(r0, v0, n, t)`
Analytically propagates a deputy spacecraft's relative state in the LVLH (RSW) rotating frame using the Clohessy-Wiltshire state transition matrix. Valid for near-circular chief orbits and small separations.

**Inputs**
- `r0` — initial relative position `[R;S;W]` (km, 3x1)
- `v0` — initial relative velocity, CW rotating-frame derivative (km/s, 3x1) — **not** inertial relative velocity (see file header for conversion)
- `n` — chief mean motion (rad/s)
- `t` — time vector (s)

**Outputs**
- `r_hist` — relative position history (km, 3xN)
- `v_hist` — relative velocity history, rotating-frame (km/s, 3xN)

### `result = cwRendezvous(dr0, dv0, n, tf, varargin)`
Computes a two-impulse CW rendezvous solution that transfers a deputy from an initial relative state to a target final relative state (default: origin/zero velocity) in transfer time `tf`.

**Inputs**
- `dr0, dv0` — initial relative position (km)/velocity (km/s, CW frame), 3x1
- `n` — chief mean motion (rad/s), `tf` — transfer time (s)
- `drf, dvf` — optional desired final relative state (default zero)
- Note: avoid `tf = k*(T/2)` (singular `Phi_rv`); prefer e.g. `0.4*T, 0.6*T, 0.75*T`

**Outputs (struct fields)**
- `dv1, dv2` (3x1 km/s), `dv1_mag, dv2_mag, total_dv_m_s`
- `r_hist, v_hist, t_hist` — 100-point trajectory through the maneuver

---

## Orbit Design

### `[result, varargout] = orbitLifetime(orb, varargin)`
Estimates orbital lifetime due to atmospheric drag using an orbit-averaged propagation (suitable for month/year timescales, not high-fidelity).

**Inputs**
- `orb` — struct from `earthOrbit()`
- Options: `'CdAm'` (0.01 m^2/kg), `'MaxYears'` (30), `'StepOrbits'` (10), `'Plot'` (false)

**Outputs (struct fields)**
- `t_days, a_km, e, alt_peri_km, alt_apo_km` — history vectors
- `lifetime_days`/`lifetime_years`, `reentered` (bool), `CdAm`
- `fig` — optional second output if `'Plot'` true

### `frz = frozenOrbit(alt_km, i_deg, varargin)`
Computes the frozen eccentricity and argument of perigee at which J3-induced forced eccentricity exactly balances J2-driven apsidal rotation (zero secular drift in `e` and `omega`).

**Inputs**
- `alt_km` — circular altitude (km), `i_deg` — inclination (deg)
- `'omega_deg'` — 90 (perigee north, default) or 270 (perigee south)

**Outputs (struct fields)**
- `a, e, i, omega, p` (km/deg), `alt_peri, alt_apo` (km), `period` (s)
- `e_dot_secular`, `RAAN_dot`, `omega_dot_J2`, `omega_dot_J2J3` (deg/day)

### `result = repeatingGroundTrack(N_rev, N_day, inc_deg, varargin)`
Iteratively solves (with J2 corrections) for the semi-major axis of an orbit that repeats its ground track every `N_rev` revolutions in `N_day` sidereal days.

**Inputs**
- `N_rev, N_day` — repeat-cycle revolutions/sidereal days, `inc_deg` — inclination (deg)
- `'e'` — eccentricity (default 0)

**Outputs (struct fields)**
- `a_km, alt_km, e, i_deg, N_rev, N_day, period_s`
- `RAAN_dot_deg_day`, `ground_track_spacing_deg`, `orb`

### `varargout = lvlhFrame(arg1, arg2, varargin)`
Transforms position/velocity vectors between ECI and the LVLH (RSW: radial/along-track/cross-track) frame, including inverse and relative (chaser-vs-target) transforms.

**Usage**
- `[r_lvlh, v_lvlh] = lvlhFrame(r_eci, v_eci)` — ECI → LVLH
- `DCM = lvlhFrame(r_eci, v_eci, 'dcm')` — ECI→LVLH rotation matrix only
- `[r_eci, v_eci] = lvlhFrame(r_lvlh, v_lvlh, r_ref_eci, v_ref_eci, 'inverse')`
- `[dr_lvlh, dv_lvlh] = lvlhFrame(r_chaser_eci, v_chaser_eci, r_target_eci, v_target_eci, 'relative')`

### `result = betaAngle(orb, duration_days, varargin)`
Computes the beta angle (angle between orbital plane and Sun direction, `beta = asin(dot(h_hat, sun_hat))`) over time, accounting for J2 RAAN precession and Earth's heliocentric motion, plus per-orbit eclipse statistics.

**Inputs**
- `orb` — struct from `earthOrbit()`, `duration_days` — analysis length (days)
- `'StepDays'` — time step (default 1.0)

**Outputs (struct fields)**
- `t_days, beta_deg, f_eclipse, t_eclipse_min, RAAN_deg` — history vectors
- `beta_max_deg, beta_min_deg, eclipse_free` (bool), `orb`

### `fig = plotBetaAngle(result)`
Two-panel plot of beta-angle history (with eclipse zone shaded) and eclipse duration per orbit vs. time.

**Inputs**
- `result` — output of `betaAngle()`

**Outputs**
- `fig` — figure handle

---

## Constellation Design

### `[sats, info] = walkerConstellation(inc_deg, alt_km, T, P, F, varargin)`
Generates a Walker delta constellation (T/P/F notation) as an array of `earthOrbit` structs distributed across `P` planes with `T/P` satellites each.

**Inputs**
- `inc_deg` — inclination (deg), `alt_km` — circular altitude (km)
- `T` — total satellites, `P` — number of planes (must divide `T`), `F` — phasing parameter (0 <= F <= P-1)
- `'epoch_jd'` — epoch for all orbits (default J2000)

**Outputs**
- `sats` — 1xT struct array of `earthOrbit` structs, each with extra fields `plane_idx, sat_idx_in_plane, sat_idx_total`
- `info` — struct of constellation parameters (e.g. satellites/plane, period)

### `cov = constellationCoverage(sats, duration_hr, varargin)`
Combined coverage analysis for an entire constellation over a lat/lon grid — analogous to `coverageAnalysis` but aggregating visibility across all satellites.

**Inputs**
- `sats` — 1xT struct array of `earthOrbit` structs
- `duration_hr` — analysis duration (hours)
- Options: `'MinElevation'` (10°), `'GridRes'` (5°), `'StepSize'` (60 s)

**Outputs (struct fields)**
- `coverage_frac, revisit_mean_hr, revisit_max_hr, n_passes` (NlatxNlon, compatible with `plotCoverage`)
- `lat_vec, lon_vec, duration_hr, orb, min_elevation_deg, n_sats`

### `fig = plotConstellationGroundTrack(sats, varargin)`
Plots ground tracks for all satellites in a constellation, optionally colored by orbital plane.

**Inputs**
- `sats` — 1xT struct array (e.g. from `walkerConstellation`)
- Options: `'NumOrbits'` (1), `'J2'` (true), `'ColorByPlane'` (true), `'GroundStations'`, `'Title'`

**Outputs**
- `fig` — figure handle

---

## Communications & Sensors

### `[result, varargout] = linkBudget(varargin)`
Two-way RF link budget analysis, either as a static snapshot at a given slant range or as a full time-series over a satellite pass (3 orbital periods).

**Usage**
- Static: `result = linkBudget(Name, Value, ...)`, optionally `[result, fig] = linkBudget(..., 'Plot', true)`
- Pass mode: `result = linkBudget(orb, gs, Name, Value, ...)` — triggered when arg1 has field `.a` (orbit) and arg2 has field `.lat` (ground station)

**Key inputs (static mode)**
- `'Freq_GHz'` (2.0), `'P_tx_dBW'` (0), `'G_tx_dBi'`/`'G_rx_dBi'` (6/45), `'T_sys_K'` (135)
- `'DataRate_bps'` (1e6), `'ReqEbN0_dB'` (10), `'L_atm_dB'`/`'L_point_dB'`, `'Range_km'` (1000)
- Pass mode also: `gs` struct (`.lat .lon .alt .name` + optional per-station link overrides), `'MinElevation'` (5°)

**Outputs (struct fields)**
- Static: `EIRP_dBW, FSPL_dB, P_rx_dBW, C_N0_dBHz, Eb_N0_dB, LinkMargin_dB, range_km, max_range_km`
- Pass: `t_s, jd, az_deg, el_deg, range_km, range_rate_km_s, doppler_Hz, Eb_N0_dB, LinkMargin_dB, passes, static_budget, gs, orb`

### `[fp, varargout] = sensorFootprint(lat_sat, lon_sat, alt_km, half_angle_deg, varargin)`
Computes (and optionally plots) the ground footprint of a conical sensor from a satellite position, supporting off-nadir pointing.

**Inputs**
- `lat_sat, lon_sat` (deg), `alt_km` — altitude above Earth (km), `half_angle_deg` — sensor half-cone angle
- Options: `'PointAz_deg'`/`'PointEl_deg'` (boresight azimuth/off-nadir angle, default 0), `'N_pts'` (360), `'Plot'` (false), `'PlotAxes'`, `'FillColor'`/`'FillAlpha'`, `'ShowSubSat'`/`'ShowCenter'`

**Outputs (struct fields)**
- `lat_boundary, lon_boundary` (Nx1) — footprint polygon
- `lat_center, lon_center` — footprint center (off-nadir)
- `earth_central_angle_deg, footprint_radius_km, footprint_area_km2, max_range_km, min_elevation_deg`
- `fig` — optional second output if `'Plot'` true
