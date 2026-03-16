# Trajectory Tools

MATLAB toolkit for preliminary mission design using the **patched-conic approximation**. Supports Earth–Moon transfers and interplanetary transfers with body-centric departure/arrival visualizations.

---

## Overview

The patched-conic method splits a multi-body trajectory into a sequence of two-body problems stitched together at sphere-of-influence (SOI) boundaries. It is fast and accurate enough for early-phase mission analysis and trade studies.

Key capabilities:
- Hohmann / direct lunar transfers with inclination and argument-of-periapsis control
- Bi-elliptic plane-change transfers to polar lunar orbits (e.g. south pole missions)
- Interplanetary transfers with pork-chop plot and best-launch-date search
- Combined speed-change + plane-change burns (Oberth-effect exploitation)
- 3D ecliptic-frame overview and body-centric departure/arrival plots

---

## File Structure

| File | Purpose |
|------|---------|
| `constants.m` | Physical constants and body parameters (μ, radius, SOI, obliquity, orbital elements) for Earth, Moon, Mars, Venus, Jupiter, and major asteroid belt bodies |
| `patchedConicTransfer.m` | Core solver — lunar and interplanetary transfer ΔV and trajectory parameters |
| `plotPatchedConic.m` | All visualization: 3D ecliptic overview, body-centric departure/arrival plots, bi-elliptic intermediate orbits |
| `lambertSolver.m` | Universal-variable Lambert solver (Stumpff functions) |
| `porkChopPlot.m` | Pork-chop contour plot over a departure/arrival date grid |
| `findBestLaunchDate.m` | Grid search for minimum-ΔV launch window |
| `julianDate.m` | Gregorian → Julian Date conversion |
| `orbitalStateCircular.m` | Circular orbit state vector (position + velocity) |
| `example_lunar_transfer.m` | Basic Earth–Moon transfer example |
| `example_lunar_south_pole.m` | South-pole mission: direct vs bi-elliptic polar capture comparison |
| `example_interplanetary_transfer.m` | Earth–Mars (or other planet) transfer example |

---

## Quick Start

Open MATLAB, add the repository folder to your path, then run any example script:

```matlab
run('example_lunar_south_pole.m')
run('example_lunar_transfer.m')
run('example_interplanetary_transfer.m')
```

---

## Transfer Options

`patchedConicTransfer(departBody, arrivalBody, options)` accepts:

| Option | Default | Description |
|--------|---------|-------------|
| `departureAltitude` | 200 km | Parking orbit altitude at departure body |
| `arrivalAltitude` | 100 km | Periapsis altitude at arrival body |
| `arrivalApogeeAltitude` | = `arrivalAltitude` | Apoapsis altitude of capture orbit (set > `arrivalAltitude` for elliptical capture) |
| `departureInclination` | 0° | Parking orbit inclination (deg); drives plane-change ΔV at TLI |
| `arrivalInclination` | 0° | Capture orbit inclination (deg) |
| `arrivalArgOfPeriapsis` | 0° | Argument of periapsis ω (deg); 90° places perilune above the north pole |
| `transferMode` | `'direct'` | `'direct'` — single combined LOI + plane change; `'biElliptic'` — three-burn sequence |
| `biEllipticApoapsisAltitude` | 0 (→ Moon SOI) | Intermediate apoapsis altitude for bi-elliptic transfers; 0 uses Moon SOI (maximum savings) |

---

## Bi-elliptic Plane Change

For high-inclination captures (especially polar orbits), a direct combined burn at perilune is ΔV-expensive because the plane change is performed at maximum speed. The bi-elliptic sequence moves the plane change to the highest—and therefore slowest—point of an intermediate orbit:

1. **Burn 1 (perilune):** Hyperbola → equatorial ellipse at high apoapsis. No plane change.
2. **Burn 2 (apolune):** Pure plane-change burn at minimum speed. ΔV = 2 · v_apo · sin(Δi / 2).
3. **Burn 3 (perilune):** Lower apoapsis from intermediate value to final capture orbit.

Savings grow monotonically with the intermediate apoapsis radius. At the Moon's SOI (~66 200 km) a 90° plane change saves roughly 0.6–1.0 km/s compared to the direct approach.

---

## Lunar South Pole Mission

`example_lunar_south_pole.m` designs a transfer to a polar capture orbit optimized for south-pole surface access:

- **Parking orbit:** 200 km LEO at 28.5° inclination (KSC latitude)
- **Capture orbit:** 100 × 5 000 km polar ellipse, ω = 90°
  - Perilune above north pole — descent approach targets the south pole
  - High apolune reduces orbital maintenance ΔV and extends south-pole communications geometry
- Compares direct and bi-elliptic transfers side-by-side with a ΔV breakdown bar chart and 3D body-centric plots showing all intermediate orbits

---

## Dependencies

- MATLAB R2019b or later (uses `deg2rad`, `rad2deg`, `surf`, `patch`, `quiver3`)
- No additional toolboxes required
