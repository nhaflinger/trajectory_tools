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
- Hyperbolic orbit propagation for interstellar/flyby objects (e > 1) via hyperbolic Kepler's equation
- No Aerospace Toolbox required — uses JPL Standish (1992) elements

**Interstellar Intercept**
- Hyperbolic ephemeris for 1I/'Oumuamua (e = 1.2, v∞ = 26.3 km/s)
- Jupiter flyby + solar Oberth maneuver trajectory design with thermal perihelion constraint
- Design space scan over launch window × Jupiter flyby date
- Oberth ΔV trade: heliocentric v∞ vs. burn size at PSP perihelion limit

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
| `tisserandGraph.m` | Tisserand parameter graph — iso-v∞ contours per body on the (r_p, r_a) plane |
| `porkChopSequence.m` | Per-leg pork-chop plots for a multi-body flyby sequence |
| `resonantOrbits.m` | Resonant return-orbit analysis: period ratios, min v∞, apsis distances |
| `lowThrustSpiral.m` | Edelbaum analytical ΔV + tangential-thrust RK4 spiral with propellant accounting |
| `plotLowThrustSpiral.m` | Spiral visualization: trajectory, altitude profile, and mass history |
| `lowThrustInterplanetary.m` | Three-phase low-thrust interplanetary budget: departure spiral + heliocentric + arrival capture |
| `plotLowThrustInterplanetary.m` | Four-panel trajectory visualization: heliocentric cruise, departure/arrival spirals, mass history |
| `porkChopLowThrust.m` | Low-thrust pork-chop: propellant fraction vs. launch date/TOF with impulsive overlay |
| `julianDate.m` | Gregorian → Julian Date conversion |
| `example_lunar_transfer.m` | Basic Earth–Moon transfer example |
| `example_lunar_south_pole.m` | South-pole mission: direct vs bi-elliptic polar capture comparison |
| `example_interplanetary_transfer.m` | Earth–Mars transfer with Lambert solver and pork-chop plot |
| `example_evj.m` | Earth–Venus–Jupiter gravity-assist trajectory with direct EJ comparison |
| `example_emej_europa_clipper.m` | Europa Clipper-style EMEJ trajectory (MEGA: Mars + Earth flybys) with direct EJ comparison |
| `example_ga_explorer.m` | Gravity-assist scenario explorer: compares six E→J flyby architectures using all three new tools |
| `example_low_thrust.m` | Low-thrust orbit transfer examples: LEO→GEO, LEO→lunar distance, propulsion trade study |
| `example_low_thrust_interplanetary.m` | Earth→Mars low-thrust budget (LV provides Earth escape), impulsive comparison, and pork-chop |
| `example_dawn.m` | NASA Dawn mission reproduction: Earth→Mars flyby→Vesta→Ceres with NSTAR ion engine budget |
| `example_oumuamua.m` | Interstellar intercept design: Earth→Jupiter flyby→Solar Oberth→1I/'Oumuamua |

---

## Quick Start

Open MATLAB, add the repository folder to your path, then run any example script:

```matlab
run('example_interplanetary_transfer.m')
run('example_lunar_south_pole.m')
run('example_lunar_transfer.m')
run('example_evj.m')
run('example_emej_europa_clipper.m')
run('example_ga_explorer.m')        % multi-architecture E→J comparison
run('example_low_thrust.m')         % low-thrust spiral transfers and propulsion trade study
run('example_dawn.m')               % NASA Dawn mission reproduction
run('example_oumuamua.m')           % interstellar intercept design: Jupiter flyby + solar Oberth
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
| **1I/'Oumuamua** | **Hyperbolic (e = 1.2, perihelion Sep 2017)** | — | — |

J2000 elements (Ω, ω, M₀) enable accurate eccentric, inclined orbit propagation without the Aerospace Toolbox. Hyperbolic bodies (e ≥ 1) use `t_peri_jd` (Julian Date of perihelion) instead of M₀ and propagate via the hyperbolic Kepler equation.

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

### Design workflow

The three analysis tools are meant to be used together in sequence, moving from broad topology to refined trajectory:

**1. Establish the design space — `tisserandGraph` + `resonantOrbits`**

Before running any optimization, generate the Tisserand graph for your flyby bodies to understand what gravity-assist architectures are geometrically possible. Add resonant-orbit markers to identify where repeat flybys are accessible.

```matlab
bodies = constants();
[fig, ax] = tisserandGraph({bodies.Venus, bodies.Earth, bodies.Jupiter});
resonantOrbits(bodies.Venus, 'ax', ax, 'print', false);
resonantOrbits(bodies.Earth, 'ax', ax, 'print', false);
```

Read the graph to identify candidate sequences: a good sequence is one where each flyby body's contours lie between the departure orbit and the destination, allowing the spacecraft to climb from one to the other with minimal powered burns.

**2. Scan candidate windows — `porkChopSequence`**

For each candidate sequence, generate per-leg pork-chop plots to see when good opportunities exist. Each subplot shows v∞ at arrival for one leg; compatible windows appear where the v∞ at the end of leg i matches the v∞ at the start of leg i+1.

```matlab
bSeq      = {bodies.Earth, bodies.Venus, bodies.Jupiter};
tofRanges = [100, 250; 600, 1200];
best      = findBestFlybyWindow(bSeq, jdStart, jdEnd, tofRanges);
porkChopSequence(bSeq, jdStart, jdEnd, tofRanges, struct('markBest', best));
```

The ★ marker shows the grid-search optimum. Use the plots to check whether the optimum sits in an isolated minimum-energy island (robust window) or on a broad shallow slope (sensitive to timing errors).

**3. Compute and compare — `flybySequence` + Tisserand overlay**

Compute the full trajectory at the best window and overlay it on the Tisserand graph to verify the gravity assists are doing the expected work.

```matlab
result = flybySequence(bSeq, best.departureJD, best.tofDays, seqOpts);
tisserandGraph({bodies.Venus, bodies.Earth, bodies.Jupiter}, ...
    struct('trajectory', result));
```

Check that consecutive leg-orbit points sit on the same v∞ contour of each flyby body — this confirms the gravity assist is free (no powered flyby burn required).

**4. Compare architectures — `example_ga_explorer.m`**

Run the explorer script to evaluate multiple sequences side by side and identify which architecture best suits your mission's ΔV and TOF requirements:

```matlab
run('example_ga_explorer.m')
```

---

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

---

## Tisserand Graph

`tisserandGraph(bodies)` plots iso-v∞ contours for each body in the `(r_p, r_a)` plane (periapsis vs. apoapsis distance in AU).

```matlab
bodies = constants();
bSeq   = {bodies.Earth, bodies.Mars, bodies.Jupiter};
[fig, ax] = tisserandGraph(bSeq);
```

To overlay a computed flyby sequence on the graph:

```matlab
result = flybySequence(bSeq, departJD, tofDays);
tisserandGraph(bSeq, struct('trajectory', result));
```

### Interpreting the Tisserand Graph

**The axes**

Each point `(r_p, r_a)` on the graph represents a class of heliocentric orbits — all ellipses with that periapsis and apoapsis distance, regardless of orientation. The dashed diagonal `r_p = r_a` is the circular-orbit locus; planets sit on it as dots.

**The contour lines**

For each body, the colored contour families show all heliocentric orbits that encounter that body with a given v∞. The contour values come from the Tisserand parameter:

```
T    = a_planet/a_spacecraft + 2·√((a_sc/a_planet)·(1−e²))
v∞   = v_planet · √(3 − T)
```

A contour labeled `2 km/s` for Earth means: any spacecraft on an orbit touching that line arrives at Earth with exactly 2 km/s hyperbolic excess speed. Lower-numbered contours (closer to each planet's dot) correspond to more efficient, lower-energy encounters.

**The key insight: gravity assists move you along a contour**

When the spacecraft flies by a planet with no thrust:
- The v∞ **magnitude** is conserved
- The **direction** changes, reshaping the heliocentric orbit
- On the Tisserand graph this appears as **sliding along that planet's contour** to a new `(r_p, r_a)` point

A **powered flyby** (burn at periapsis via Oberth effect) lets you jump between contours of the same body at the cost of some ΔV, applied at maximum speed.

**Reading a trajectory sequence**

Each leg's transfer ellipse maps to one point `(r_p, r_a)`. The trajectory overlaid by `example_ga_explorer.m` appears as a dashed colored path through the graph:

```
Leg 1 point  →  flyby at planet X  →  Leg 2 point  →  flyby at planet Y  →  ...
```

At each flyby, the jump from one leg's point to the next **should lie along that planet's v∞ contour**. If both points fall on the same contour value → the gravity assist is free. If they fall on different contours → a powered flyby burn was required.

**Reaching a distant destination efficiently** means climbing the `r_a` axis with as little ΔV as possible. A well-designed gravity-assist sequence looks like:

1. Depart the inner planet with modest v∞ (dot near departure body's line, modest `r_a`)
2. Each flyby slides the dot along a contour, raising `r_a` toward the destination
3. The final dot lands on a low v∞ contour of the destination → small arrival burn

**What to look for**

| You observe... | It means... |
|----------------|-------------|
| Two consecutive trajectory dots on the same body's contour | Free gravity assist — no ΔV spent |
| Dots on widely separated contours of the flyby body | Powered flyby needed, or poor sequence choice |
| Path climbing smoothly up the `r_a` axis | Each assist efficiently raises the apoapsis toward the target |
| Final dot on a low-numbered destination contour | Low arrival v∞ → small capture burn |
| × marker near a trajectory dot | A repeat flyby at that body is geometrically accessible |

**Resonant × markers**

The `×` symbols placed by `resonantOrbits` mark specific `(r_p, r_a)` points where the spacecraft's orbital period is a simple integer ratio p:q of the flyby planet's period. After the flyby the spacecraft completes p revolutions, the planet completes q, and they meet again at the same point — enabling a second flyby with no extra targeting. These are the anchor points for v∞ leveraging (VILT) strategies.

### `tisserandGraph` options

| Option | Default | Description |
|--------|---------|-------------|
| `vInfValues` | `[0.5 1 2 3 5 10 20]` | v∞ contour levels (km/s) |
| `nGrid` | `400` | Grid resolution |
| `rpLim` | auto | Periapsis range `[min max]` AU |
| `raLim` | auto | Apoapsis range `[min max]` AU |
| `trajectory` | — | `flybySequence` result to overlay |
| `showLabels` | `true` | Label contour v∞ values |

### Legend guide

| Symbol | Meaning |
|--------|---------|
| ● colored dot | Planet location on its circular orbit |
| — colored line (labeled) | Iso-v∞ contour — the number is v∞ in km/s |
| × colored mark | p:q resonant return orbit — label gives the ratio and minimum v∞ |
| white dashed path | Trajectory sequence connecting leg-orbit points |
| ○ numbered dot (L1, L2, …) | Transfer-leg orbit position in (r_p, r_a) space |
| diagonal dashed line | Circular orbit locus (r_p = r_a) — planets sit on this line |

**Reading a flyby:** at each intermediate body, the jump from one L# dot to the next should lie on one of that body's contours at the same labeled v∞ value.  If both dots are on the same contour → free gravity assist.  If they are on different contours → powered flyby ΔV was required.

### Other Tisserand graph variants

The `(r_p, r_a)` graph implemented here is the most common form for trajectory topology, but several other projections are used in mission design, each optimized for a different question:

**v∞ vs. orbital period (period-Tisserand)**

X-axis: v∞ at the flyby body (km/s); Y-axis: spacecraft orbital period (years).  As v∞ increases, a wider range of periods becomes accessible from a flyby of that body.  The key advantage of this projection is that **resonant orbits appear as horizontal lines** at exact period ratios — 1 yr, 2 yr, 1.5 yr, etc. — making v∞ leveraging (VILT) chains immediately visible.  Used extensively in designing repeated flyby sequences around the outer planets.

**v∞ vs. semi-major axis (energy Tisserand)**

X-axis: v∞ at the flyby body (km/s); Y-axis: spacecraft semi-major axis (AU).  Equivalent to the period version via Kepler's third law.  Shows directly how the accessible semi-major axis range expands with increasing v∞, useful for bounding the ΔV required to reach a target orbit from a given flyby speed.

**Linked v∞ graph (bi-body or Tisserand-Poincaré graph)**

X-axis: v∞ at body 1 (km/s); Y-axis: v∞ at body 2 (km/s).  For a trajectory that encounters both bodies, each heliocentric transfer ellipse maps to a curve in this space.  A flyby at body 1 moves the spacecraft along a **vertical line** (v∞₁ conserved, orbit reshapes); a flyby at body 2 moves it along a **horizontal line** (v∞₂ conserved).  The intersection of curves from two different bodies shows which (v∞₁, v∞₂) pairs are achievable from the same ellipse — directly identifying free-gravity-assist connections.  This is the standard graph for designing multi-body VILT chains (Campagnola & Russell 2010).

**Choosing a projection**

| Question | Best graph |
|----------|-----------|
| What flyby sequences are geometrically possible? | r_p vs. r_a (this implementation) |
| Where are the resonant return windows? | v∞ vs. period (period-Tisserand) |
| Can I connect two bodies with a free gravity assist? | v∞₁ vs. v∞₂ (bi-body / linked graph) |
| How much does a flyby change my orbital energy? | v∞ vs. semi-major axis (energy-Tisserand) |

---

## Resonant Orbit Analysis

`resonantOrbits(body)` enumerates p:q resonant return orbits for a flyby planet — the set of heliocentric orbits that return the spacecraft to the planet after exactly p spacecraft revolutions per q planet revolutions.

```matlab
bodies = constants();
resonantOrbits(bodies.Earth);               % print table
resonantOrbits(bodies.Mars);

% Overlay on a Tisserand graph
fig = tisserandGraph({bodies.Earth, bodies.Mars, bodies.Jupiter});
resonantOrbits(bodies.Earth, 'ax', gca);
resonantOrbits(bodies.Mars,  'ax', gca);
```

**Physics:**  For resonance p:q, the spacecraft semi-major axis is `a_sc = a_p·(q/p)^(2/3)`.  The minimum v∞ at a tangent (minimum-energy) encounter is:

```
v∞_min = v_p · |√(2 − a_p/a_sc) − 1|
```

Outward resonances (`a_sc > a_p`) have periapsis at the planet; inward resonances (`a_sc < a_p`) have apoapsis at the planet.

### `resonantOrbits` options

| Option / name-value | Default | Description |
|---------------------|---------|-------------|
| `maxN` | `5` | Maximum integer in p and q |
| `ax` | — | Axes handle — marks resonances on Tisserand graph |
| `print` | `true` | Print summary table |

---

## Per-Leg Pork-Chop Sequence

`porkChopSequence(bodies, jdStart, jdEnd, tofRanges)` creates one pork-chop subplot per leg showing v∞ at arrival as a function of departure date and TOF.  The departure window for leg i is automatically offset from `jdStart/jdEnd` by the accumulated minimum TOF of the preceding legs.

```matlab
bodies = constants();
bSeq   = {bodies.Earth, bodies.Venus, bodies.Jupiter};

jdStart   = julianDate(2026, 1, 1);
jdEnd     = julianDate(2030, 12, 31);
tofRanges = [100, 250;    % leg 1: E→V
             600, 1200];  % leg 2: V→J

porkChopSequence(bSeq, jdStart, jdEnd, tofRanges);
```

To mark the optimum from `findBestFlybyWindow`:

```matlab
best = findBestFlybyWindow(bSeq, jdStart, jdEnd, tofRanges);
porkChopSequence(bSeq, jdStart, jdEnd, tofRanges, struct('markBest', best));
```

### `porkChopSequence` options

| Option | Default | Description |
|--------|---------|-------------|
| `nDep` | `60` | Departure-date grid points per leg |
| `nTof` | `60` | TOF grid points per leg |
| `transferTypes` | `'type1'` each | Per-leg `'type1'` or `'type2'` |
| `markBest` | — | Struct `.departureJD` + `.tofDays` to mark ★ |
| `cLimPct` | `[5 85]` | Colorscale percentile clamp |

---

### Flyby physics

At each intermediate body the incoming and outgoing v∞ vectors come directly from the adjacent Lambert solutions.

**Free gravity assist** (`|v∞_in| = |v∞_out|`):

sin(δ/2) = 1 / (1 + r_p · v∞² / μ) → r_p = (μ/v∞²) · (1/sin(δ/2) − 1)

**Powered flyby** (`|v∞_in| ≠ |v∞_out|`) — ΔV applied at periapsis via Oberth effect:

ΔV = |√(v∞_out² + 2μ/r_p) − √(v∞_in² + 2μ/r_p)|

Savings from gravity assists grow with the flyby body's mass and with higher incoming v∞.  For inner-planet assists (Venus, Earth) the deflection angle is large but the speed boost is modest; for Jupiter flybys the speed change can exceed several km/s.

---

## Low-Thrust Spiral Transfers

`lowThrustSpiral(body, r0_km, r1_km)` computes the ΔV and propellant budget for a continuous low-thrust transfer between two circular orbits.

```matlab
bodies = constants();

% LEO -> GEO with a Hall thruster (no plane change)
opts           = struct();
opts.thrustN   = 0.300;   % N
opts.isp       = 1800;    % s (Hall thruster)
opts.wetMass   = 500;     % kg

res = lowThrustSpiral(bodies.Earth, ...
    bodies.Earth.radius + 200,    ...   % 200 km LEO
    bodies.Earth.radius + 35786,  ...   % GEO
    opts);

fprintf('dV = %.3f km/s,  TOF = %.0f days,  Prop = %.0f kg\n', ...
    res.deltaV, res.tofDays, res.propellantMass);

plotLowThrustSpiral(res, bodies.Earth);
```

### Physics

**Edelbaum analytical ΔV** (1961): optimal continuous-thrust spiral between two circular orbits with simultaneous inclination change:

```
ΔV = √(v₀² + v₁² − 2·v₀·v₁·cos(π·Δi/2))
```

- `v₀`, `v₁` — circular orbital speeds at `r0`, `r1`
- `Δi` — inclination change in **radians**
- For `Δi = 0`: reduces to `|v₀ − v₁|` (correct result for a slow spiral)
- For combined orbit-raising + plane change, Edelbaum distributes the plane change optimally throughout the spiral, giving **lower ΔV than a two-burn impulsive sequence** when `Δi` is large

**Oberth-effect trade-off**: low-thrust spirals require *more* total ΔV than equivalent impulsive burns (the Oberth advantage of burning at maximum speed is lost), but the high Isp of electric thrusters results in far *less propellant mass* consumed.

| Mission | Impulsive ΔV | Edelbaum ΔV | Note |
|---------|-------------|-------------|------|
| LEO→GEO, Δi=0° | 3.9 km/s | 4.7 km/s | Spiral costs +0.8 km/s |
| LEO→GEO, Δi=28.5° | ~6.3 km/s | ~6.0 km/s | Spiral is *cheaper* — plane change advantage |

**RK4 tangential-thrust propagator**: for `thrustN > 0`, the actual spiral trajectory is propagated in 2-D Cartesian coordinates with thrust always along the velocity vector (prograde for orbit-raising, retrograde for deorbit). This gives accurate trajectory geometry at low thrust-to-weight ratios.

**Propellant accounting** via Tsiolkovsky: `Δm = m₀(1 − exp(−ΔV/(Isp·g₀)))`.

### `lowThrustSpiral` options

| Option | Default | Description |
|--------|---------|-------------|
| `inclinationChangeDeg` | `0` | Combined inclination change (deg) — Edelbaum handles optimally |
| `thrustN` | `0` | Engine thrust (N); `0` computes Edelbaum ΔV only (no trajectory) |
| `isp` | `3000` | Specific impulse (s) |
| `wetMass` | `1000` | Initial spacecraft wet mass (kg) |
| `nStepsPerOrbit` | `50` | RK4 steps per osculating orbital period |
| `nOutputPoints` | `2000` | Decimated trajectory output points |

### Output fields

| Field | Description |
|-------|-------------|
| `result.deltaV` | Edelbaum ΔV (km/s) |
| `result.tof` | Estimated time of flight (s) — requires `thrustN > 0` |
| `result.tofDays` | Same in days |
| `result.propellantMass` | Propellant consumed (kg) |
| `result.finalMass` | Final spacecraft dry+residual mass (kg) |
| `result.details` | Component values: `v0`, `v1`, `r0`, `r1`, `thrust`, `isp`, `wetMass` |
| `result.trajectory` | Decimated RK4 history: `.t` `.x` `.y` `.r` `.speed` `.mass` |

---

## Low-Thrust Interplanetary Transfers

### Three-phase mission model

`lowThrustInterplanetary(departBody, arrivalBody, departJD, tofDays)` computes the complete ΔV and propellant budget for a low-thrust interplanetary mission as three sequential Edelbaum spirals:

| Phase | Transfer | Central body |
|-------|----------|-------------|
| 1 | Parking orbit → departure body SOI | Departure planet |
| 2 | Heliocentric cruise (dep orbit → arr orbit) | Sun |
| 3 | Arrival body SOI → parking orbit | Arrival planet |

Mass is propagated sequentially — each phase starts with the mass remaining after the previous one.

```matlab
bodies = constants();
opts = struct('thrustN', 0.1, 'isp', 3000, 'wetMass', 1000, ...
              'departureAltitude', 200, 'arrivalAltitude', 400);

res = lowThrustInterplanetary(bodies.Earth, bodies.Mars, ...
          julianDate(2026,1,1), 800, opts);

fprintf('Total ΔV: %.2f km/s   Propellant: %.0f kg   TOF: %.0f days\n', ...
    res.deltaV, res.propellantMass, res.tofDays);
```

Actual planet positions (eccentric orbits) are used for the heliocentric phase, so orbital eccentricity (e.g. Mars e = 0.093) contributes to the result. The Edelbaum formula is a lower bound — it is exact for infinitely slow spirals and becomes optimistic at higher thrust-to-weight ratios where gravity losses appear.

**The Oberth trade-off**: Low-thrust requires significantly more total ΔV than impulsive (the Oberth effect is lost), but the high Isp of electric thrusters means far less propellant mass:

| Metric | Impulsive (Isp 450 s) | Low-thrust (Isp 3000 s) |
|--------|----------------------|------------------------|
| Earth→Mars ΔV | ~4.5 km/s | ~15 km/s (all three spirals) |
| Propellant fraction | ~65% | ~40% |
| TOF | ~300 days | ~500–900 days |

### Low-thrust pork-chop

`porkChopLowThrust(departBody, arrivalBody, departJD, tofDays)` sweeps a (departure date × TOF) grid and plots the propellant mass fraction for the complete three-phase mission. White dashed contours overlay the equivalent impulsive Lambert ΔV for comparison.

```matlab
depDates = linspace(julianDate(2026,1,1), julianDate(2028,12,31), 60);
tofRange = linspace(400, 1400, 50);

porkChopLowThrust(bodies.Earth, bodies.Mars, depDates, tofRange, ...
    struct('isp', 3000, 'showImpulsive', true));
```

The launch-window structure in the low-thrust pork-chop is driven by planetary eccentricity (orbital speed varies around the orbit) and is much flatter than the impulsive equivalent — low-thrust missions are less sensitive to launch date but more sensitive to total TOF budget.

### `porkChopLowThrust` options

| Option | Default | Description |
|--------|---------|-------------|
| `isp` | `3000` | Specific impulse for propellant calculation (s) |
| `wetMass` | `1000` | Spacecraft wet mass (kg) |
| `departureAltitude` | `200` | Departure parking orbit altitude (km) |
| `arrivalAltitude` | `400` | Arrival parking orbit altitude (km) |
| `colorMode` | `'propFrac'` | `'propFrac'` propellant fraction; `'deltaV'` total ΔV |
| `cLimPct` | `[5 85]` | Colorscale percentile clamp |
| `showImpulsive` | `true` | Overlay impulsive Lambert ΔV contours (white dashed) |
| `impulsiveContours` | auto | ΔV levels (km/s) for the impulsive overlay |

---

## Interstellar Intercept: 1I/'Oumuamua

`example_oumuamua.m` designs an imaginary pre-discovery mission to intercept the first known interstellar object, 1I/'Oumuamua, using a Jupiter gravity assist and solar Oberth maneuver.

### Mission architecture

```
Earth departure (high C3)
    → Jupiter gravity assist   [bend trajectory sunward, boost heliocentric speed]
        → Solar perihelion     [Oberth kick — chemical burn at closest Sun approach]
            → Heliocentric escape toward Oumuamua
```

### Thermal constraint — Parker Solar Probe perihelion limit

The spacecraft's minimum solar distance is constrained to no closer than the Parker Solar Probe's planned final perihelion (9.86 R☉ ≈ 0.046 AU). This sets the Oberth perihelion distance and therefore the maximum achievable speed boost.

### Why the Oberth effect is so powerful here

At 0.046 AU the Sun's escape velocity is ~196 km/s. A spacecraft arriving from Jupiter on a diving ellipse reaches perihelion at ~195 km/s — nearly equal to escape velocity, so even a small kick ΔV produces a large increase in heliocentric v∞:

```
v∞² = (v_peri + ΔV)² − v_esc²
```

| Oberth ΔV | Post-burn heliocentric v∞ |
|-----------|--------------------------|
| ~0 km/s   | ~0 km/s (bound) |
| 2 km/s    | ~20 km/s |
| 5 km/s    | ~45 km/s |
| 10 km/s   | ~70 km/s |

Oumuamua's own heliocentric v∞ is 26.3 km/s — achievable with an Oberth burn of roughly 3–4 km/s at the PSP perihelion limit.

### Design space output (3 figures)

**Figure 1 — Design space scan**: three side-by-side contour plots over the (launch date 2011–2016) × (Jupiter flyby date 2012–2017) grid:
- Launch C₃ (km²/s²)
- Post-flyby perihelion (AU) with PSP thermal limit shown as an orange dashed contour
- Post-Oberth heliocentric v∞ (km/s) with Oumuamua's 26.3 km/s shown as a green contour

The white ★ marks the best trajectory found by the grid search (maximize v∞ subject to alignment and C₃ constraints).

**Figure 2 — Oberth trade curve**: heliocentric v∞ vs. Oberth ΔV at the PSP perihelion for the best trajectory. Shows the minimum ΔV to escape the Sun and the ΔV required to match Oumuamua's speed.

**Figure 3 — Solar system view**: ecliptic top-down showing:
- Oumuamua's hyperbolic path (inbound dimmer, outbound brighter, perihelion marked)
- Spacecraft: Earth→Jupiter Lambert arc, Jupiter→perihelion dashed arc, post-Oberth escape arrow
- PSP perihelion constraint circle, Sun, planet orbit traces, key date labels

### Hyperbolic ephemeris extension

`orbitalState.m` was extended to support hyperbolic bodies (e ≥ 1). Add any hyperbolic flyby object to `constants.m` with these fields:

```matlab
body.a          = -q / (e - 1);   % km, negative by convention
body.e          = 1.xxxxx;        % eccentricity > 1
body.inclination = ...;           % deg, ecliptic J2000
body.Omega       = ...;           % deg
body.omega_peri  = ...;           % deg
body.t_peri_jd   = ...;           % Julian Date of perihelion passage
```

The propagator solves `M_h = e·sinh(H) − H` via Newton-Raphson (Battin initial guess), then converts to true anomaly and rotates to heliocentric ecliptic using the same direction cosine matrix as elliptic bodies.

---

## Dependencies

- MATLAB R2019b or later (`deg2rad`, `rad2deg`, `surf`, `patch`, `quiver3`, `contourf`)
- No additional toolboxes required
