# Trajectory Tools

MATLAB toolkit for preliminary mission design using the **patched-conic approximation**. Supports Earth–Moon transfers, interplanetary transfers with Lambert solver, pork-chop launch-window analysis, and 3D ecliptic-frame visualizations — all without any additional MATLAB toolboxes.

---

## Overview

The patched-conic method splits a multi-body trajectory into a sequence of two-body problems stitched together at sphere-of-influence (SOI) boundaries. It is fast and accurate enough for early-phase mission analysis and trade studies.

Key capabilities:

**Interplanetary**
- Lambert solver for point-to-point transfers using true eccentric, inclined planet orbits
- Pork-chop plot with best-launch-date search over a multi-year window
- Type I (short-way, < 180°) and Type II (long-way, > 180°) transfer comparison
- C3, departure/arrival v∞, ΔV budget breakdown, and TCM reserve
- 3D heliocentric overview with true elliptical planet orbits and actual Lambert transfer arc

**Gravity Assists**
- Multi-body flyby sequences with any number of intermediate flyby bodies
- Free gravity assists (|v∞_in| = |v∞_out|) and powered flybys (Oberth-effect ΔV at periapsis)
- Per-flyby: deflection angle, required periapsis radius, feasibility check, max achievable deflection
- Grid search for optimal departure date and per-leg TOFs
- 3D heliocentric visualization of all legs with flyby annotations

**Lunar**
- Hohmann / direct transfers with inclination and argument-of-periapsis control
- Bi-elliptic plane-change transfers to polar lunar orbits (south pole missions)
- Combined speed-change + plane-change burns (Oberth-effect exploitation)
- 3D body-centric departure and arrival plots with intermediate orbit visualization

**Ephemeris**
- Full Keplerian propagation from J2000 elements (Ω, ω, M₀) for all major bodies
- No Aerospace Toolbox required — uses JPL Standish (1992) elements

---

## File Structure

| File | Purpose |
|------|---------|
| `constants.m` | Physical and orbital constants for all supported bodies |
| `orbitalState.m` | Heliocentric ecliptic state vector from J2000 Keplerian elements |
| `orbitalStateCircular.m` | Wrapper — delegates to `orbitalState` when full elements are present |
| `lambertSolver.m` | Universal-variable Lambert solver (Bate-Mueller-White, Stumpff functions) |
| `patchedConicTransfer.m` | Core solver — lunar and interplanetary ΔV and trajectory parameters |
| `plotPatchedConic.m` | All visualization: 3D heliocentric overview, body-centric plots, bi-elliptic orbits |
| `porkChopPlot.m` | Pork-chop contour plot over a departure/arrival date grid |
| `findBestLaunchDate.m` | Grid search for minimum-ΔV two-body launch window |
| `gravityAssist.m` | Flyby geometry analysis: deflection angle, periapsis, powered-flyby ΔV |
| `flybySequence.m` | Multi-leg Lambert chaining with gravity-assist analysis at each flyby body |
| `findBestFlybyWindow.m` | Grid search for minimum-ΔV gravity-assist launch window |
| `plotFlybySequence.m` | 3D heliocentric visualization of multi-leg gravity-assist trajectories |
| `julianDate.m` | Gregorian → Julian Date conversion |
| `example_lunar_transfer.m` | Basic Earth–Moon transfer example |
| `example_lunar_south_pole.m` | South-pole mission: direct vs bi-elliptic polar capture comparison |
| `example_interplanetary_transfer.m` | Earth–Mars transfer with Lambert solver and pork-chop plot |
| `example_gravity_assist.m` | Earth–Venus–Jupiter gravity-assist trajectory with direct comparison |

---

## Quick Start

Open MATLAB, add the repository folder to your path, then run any example script:

```matlab
run('example_interplanetary_transfer.m')
run('example_lunar_south_pole.m')
run('example_lunar_transfer.m')
run('example_gravity_assist.m')
```

---

## Supported Bodies

| Body | Orbit Elements | μ | SOI |
|------|---------------|---|-----|
| Earth | J2000 (eccentric + inclined) | 398 600 km³/s² | — |
| Moon | circular (geocentric) | 4 903 km³/s² | ~66 200 km |
| Mercury | J2000 (e = 0.206, i = 7.0°) | 22 032 km³/s² | — |
| Venus | J2000 | 324 859 km³/s² | — |
| Mars | J2000 (e = 0.093, i = 1.85°) | 42 828 km³/s² | ~577 000 km |
| Jupiter | J2000 | 1.267 × 10⁸ km³/s² | — |
| Saturn | J2000 | 3.793 × 10⁷ km³/s² | — |
| Uranus | J2000 | 5.794 × 10⁶ km³/s² | — |
| Neptune | J2000 | 6.837 × 10⁶ km³/s² | — |
| Pluto | J2000 (e = 0.249, i = 17.1°) | 870 km³/s² | — |
| Ceres, Vesta, Pallas, Hygiea | J2000 | — | — |
| Eris, Makemake, Haumea | J2000 (approx.) | — | — |

J2000 elements (Ω, ω, M₀) enable accurate eccentric, inclined orbit propagation without the Aerospace Toolbox.

---

## Interplanetary Transfer Options

`patchedConicTransfer(departBody, arrivalBody, options)` accepts:

| Option | Default | Description |
|--------|---------|-------------|
| `departureAltitude` | 200 km | Parking orbit altitude at departure body |
| `arrivalAltitude` | 400 km | Periapsis altitude of capture orbit |
| `arrivalApogeeAltitude` | = `arrivalAltitude` | Apoapsis altitude; set higher for elliptical capture |
| `departureInclination` | 0° | Parking orbit inclination — drives plane-change ΔV at departure |
| `arrivalInclination` | 0° | Capture orbit inclination |
| `departureJD` | *(none)* | Julian Date of departure; enables Lambert solver |
| `tofDays` | *(none)* | Time of flight in days; enables Lambert solver |
| `transferType` | `'type1'` | `'type1'` short-way (< 180°); `'type2'` long-way (> 180°) |

When `departureJD` and `tofDays` are provided the solver uses true eccentric planet positions and the Lambert solution. When omitted it falls back to a Hohmann approximation.

### Output fields (interplanetary)

| Field | Description |
|-------|-------------|
| `result.deltaV` | Total ΔV including TCM reserve (km/s) |
| `result.deltaVBurns` | Departure + arrival burns only (km/s) |
| `result.details.dvDeparture` | Departure burn ΔV (km/s) |
| `result.details.dvArrival` | Arrival/capture burn ΔV (km/s) |
| `result.details.dvTCM` | TCM reserve — 2% of v∞dep, 10–50 m/s range (km/s) |
| `result.details.C3` | Characteristic energy C3 = v∞dep² (km²/s²) |
| `result.details.vInfDepart` | Departure hyperbolic excess speed (km/s) |
| `result.details.vInfArrive` | Arrival hyperbolic excess speed (km/s) |

---

## Lunar Transfer Options

`patchedConicTransfer(Earth, Moon, options)` accepts:

| Option | Default | Description |
|--------|---------|-------------|
| `departureAltitude` | 200 km | LEO parking orbit altitude |
| `arrivalAltitude` | 200 km | Perilune altitude of capture orbit |
| `arrivalApogeeAltitude` | = `arrivalAltitude` | Apolune altitude; set higher for elliptical capture |
| `departureInclination` | 0° | Parking orbit inclination relative to transfer plane — drives plane-change ΔV at TLI |
| `arrivalInclination` | 0° | Capture orbit inclination relative to lunar equator |
| `arrivalArgOfPeriapsis` | 0° | Argument of periapsis of capture orbit (90° = perilune above north pole) |
| `transferMode` | `'direct'` | `'direct'` combined LOI + plane change; `'biElliptic'` three-burn sequence |
| `biEllipticApoapsisAltitude` | 0 (Moon SOI) | Intermediate apoapsis altitude for bi-elliptic mode; 0 uses Moon's SOI for maximum savings |

### Output fields (lunar — direct)

| Field | Description |
|-------|-------------|
| `result.deltaV` | Total ΔV: TLI + capture (km/s) |
| `result.tof` | Transfer time of flight — half the ellipse period (s) |
| `result.details.dvTLI` | Trans-lunar injection ΔV including departure plane change (km/s) |
| `result.details.dvCapture` | Lunar orbit insertion ΔV including arrival plane change (km/s) |
| `result.details.vInf` | Hyperbolic excess speed at Moon SOI (km/s) |
| `result.details.transferSemiMajor` | TLI transfer ellipse semi-major axis (km) |
| `result.details.r0` | Departure parking orbit radius (km) |

### Output fields (lunar — bi-elliptic, additional fields)

| Field | Description |
|-------|-------------|
| `result.details.dvLOI` | Burn 1 — hyperbola → equatorial ellipse at perilune, no plane change (km/s) |
| `result.details.dvPlaneChange` | Burn 2 — pure plane change at intermediate apolune (km/s) |
| `result.details.dvApoapsisTrim` | Burn 3 — trim apoapsis from intermediate to final value (km/s) |
| `result.details.rBiEllipticApoapsis` | Intermediate apoapsis radius used for plane change (km) |

---

## Lambert Solver

`lambertSolver(r1, r2, tof, mu)` implements the universal-variable method from Bate, Mueller & White (1971, §5.3):

- Parameter `z` spans hyperbolic (z < 0) through elliptic (0 < z < (2π)²) solutions
- `T(z)` is monotonically decreasing — a 200-point scan locates the bracket, bisection refines it to 1 × 10⁻¹⁰ relative error
- Velocities from exact Lagrange coefficients: `f = 1 − y/R1`, `g = A√(y/μ)`, `ġ = 1 − y/R2`
- Supports short-way (`isLongWay = false`) and long-way (`isLongWay = true`) transfers

---

## Pork-Chop Plot

```matlab
bodies  = constants();
jdStart = julianDate(2026, 1, 1);
jdEnd   = julianDate(2032, 12, 31);
tofDays = linspace(120, 350, 80)';

best = findBestLaunchDate(bodies.Earth, bodies.Mars, jdStart, jdEnd, tofDays);
porkChopPlot(bodies.Earth, bodies.Mars, best.departJD, tofDays);
```

The plot uses dual calendar ticks — year boundaries (major) and month starts (minor) — with the colorscale clamped to the 5th–85th percentile of valid ΔV values so the minimum-energy island is clearly visible.

---

## Bi-elliptic Plane Change (Lunar)

For high-inclination captures, a direct combined burn at perilune is ΔV-expensive because the plane change is performed at maximum speed. The bi-elliptic sequence moves the plane change to the highest — and therefore slowest — point of an intermediate orbit:

1. **Burn 1 (perilune):** Hyperbola → equatorial ellipse at high apoapsis. No plane change.
2. **Burn 2 (apolune):** Pure plane-change burn at minimum speed. ΔV = 2 · v_apo · sin(Δi / 2).
3. **Burn 3 (perilune):** Lower apoapsis from intermediate value to final capture orbit.

Savings grow monotonically with the intermediate apoapsis. At the Moon's SOI (~66 200 km) a 90° plane change saves roughly 0.6–1.0 km/s compared to the direct approach.

---

## Lunar South Pole Mission

`example_lunar_south_pole.m` designs a transfer to a polar capture orbit optimized for south-pole surface access:

- **Parking orbit:** 200 km LEO at 28.5° inclination (KSC latitude)
- **Capture orbit:** 100 × 5 000 km polar ellipse, ω = 90°
  - Perilune above north pole — descent approach targets the south pole
  - High apolune reduces orbital maintenance ΔV and extends south-pole communications geometry
- Reports transfer orbit parameters (a, e, period), v∞ at Moon SOI, and approach speed at perilune for both modes
- TCM budget (2% of TLI ΔV, minimum 10 m/s) included in both the comparison table and the bar chart
- Compares direct and bi-elliptic transfers with a full ΔV breakdown and 3D plots showing all intermediate orbits

---

---

## Gravity Assist Design

### Quick start

```matlab
bodies = constants();
bSeq   = {bodies.Earth, bodies.Venus, bodies.Jupiter};

% Search for best departure window
tofRanges = [100, 250; 600, 1200];   % [min, max] days per leg
best = findBestFlybyWindow(bSeq, julianDate(2026,1,1), julianDate(2030,12,31), tofRanges);

% Full trajectory analysis
result = flybySequence(bSeq, best.departureJD, best.tofDays);
plotFlybySequence(result, bSeq);
```

### `flybySequence` options

| Option | Default | Description |
|--------|---------|-------------|
| `departureAltitude` | 200 km | Parking orbit altitude at departure body |
| `arrivalAltitude` | 400 km | Capture orbit periapsis altitude at arrival body |
| `arrivalApogeeAltitude` | = `arrivalAltitude` | Capture orbit apoapsis altitude |
| `departureInclination` | 0° | Departure parking orbit inclination |
| `arrivalInclination` | 0° | Capture orbit inclination |
| `flybyAltitudes` | 300 km each | 1×(N−2) periapsis altitude at each flyby body |
| `atmosphereAltitudes` | 0 km each | 1×(N−2) minimum safe altitude above body surface |
| `transferTypes` | `'type1'` each | Per-leg `'type1'` (short-way) or `'type2'` (long-way) |

### Output fields (`flybySequence`)

| Field | Description |
|-------|-------------|
| `result.deltaV` | Total ΔV: departure + powered flybys + arrival + TCM (km/s) |
| `result.deltaVBurns` | Same without TCM reserve (km/s) |
| `result.tof` | Total time of flight (days) |
| `result.details.dvDeparture` | Departure burn ΔV (km/s) |
| `result.details.dvPoweredFlybys` | Sum of powered-flyby burns (km/s); 0 for free GAs |
| `result.details.dvArrival` | Arrival/capture burn ΔV (km/s) |
| `result.details.dvTCM` | TCM reserve (km/s) |
| `result.details.C3` | Departure characteristic energy (km²/s²) |
| `result.legs(i)` | Per-leg Lambert solution (r, v vectors, v∞ in/out) |
| `result.flybys(j)` | Per-flyby geometry (deflection, periapsis, feasibility, dvPowered) |

### Flyby physics

At each intermediate body the incoming and outgoing v∞ vectors come directly from the adjacent Lambert solutions.

**Free gravity assist** (`|v∞_in| = |v∞_out|`):

sin(δ/2) = 1 / (1 + r_p · v∞² / μ) → r_p = (μ/v∞²) · (1/sin(δ/2) − 1)

**Powered flyby** (`|v∞_in| ≠ |v∞_out|`) — ΔV applied at periapsis via Oberth effect:

ΔV = |√(v∞_out² + 2μ/r_p) − √(v∞_in² + 2μ/r_p)|

Savings from gravity assists grow with the flyby body's mass and with higher incoming v∞.  For inner-planet assists (Venus, Earth) the deflection angle is large but the speed boost is modest; for Jupiter flybys the speed change can exceed several km/s.

---

## Dependencies

- MATLAB R2019b or later (`deg2rad`, `rad2deg`, `surf`, `patch`, `quiver3`, `contourf`)
- No additional toolboxes required
