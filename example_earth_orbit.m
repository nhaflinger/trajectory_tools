%% example_earth_orbit.m
% Demonstrates earthOrbit(), plotGroundTrack(), and plotOrbit3D() for
% several representative Earth orbit types.
%
% Functions used: earthOrbit, coe2eci, eci2coe, plotGroundTrack, plotOrbit3D

SCRIPT_DIR = fileparts(mfilename('fullpath'));
OUT_DIR    = fullfile(SCRIPT_DIR, 'output_earth_orbit');
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

%% ── 1. ISS-like LEO (circular, 51.6 deg inclination) ────────────────────────
iss = earthOrbit('circular', 410, 51.6);
fprintf('ISS-like LEO\n');
fprintf('  a        = %.1f km\n', iss.a);
fprintf('  alt      = %.0f km\n', iss.alt_peri);
fprintf('  i        = %.1f deg\n', iss.i);
fprintf('  period   = %.1f min\n\n', iss.period/60);

fig1 = plotGroundTrack(iss, 'NumOrbits', 3, 'J2', true);
saveas(fig1, fullfile(OUT_DIR, 'iss_groundtrack.png'));

fig2 = plotOrbit3D(iss);
saveas(fig2, fullfile(OUT_DIR, 'iss_orbit3d.png'));

%% ── 2. Sun-synchronous orbit (Earth observation) ────────────────────────────
sso = earthOrbit('sso', 550);
fprintf('Sun-Synchronous Orbit  (550 km)\n');
fprintf('  i (SSO)  = %.4f deg\n', sso.i);
fprintf('  period   = %.2f min\n\n', sso.period/60);

fig3 = plotGroundTrack(sso, 'NumOrbits', 5, 'J2', true);
saveas(fig3, fullfile(OUT_DIR, 'sso_groundtrack.png'));

fig4 = plotOrbit3D(sso);
saveas(fig4, fullfile(OUT_DIR, 'sso_orbit3d.png'));

%% ── 3. Molniya orbit ─────────────────────────────────────────────────────────
mol = earthOrbit('molniya', 0);
fprintf('Molniya Orbit\n');
fprintf('  a        = %.0f km\n', mol.a);
fprintf('  alt_apo  = %.0f km\n', mol.alt_apo);
fprintf('  alt_peri = %.0f km\n', mol.alt_peri);
fprintf('  i        = %.3f deg  (critical, no apsidal drift)\n', mol.i);
fprintf('  period   = %.1f h\n\n', mol.period/3600);

fig5 = plotGroundTrack(mol, 'NumOrbits', 4, 'J2', true);
saveas(fig5, fullfile(OUT_DIR, 'molniya_groundtrack.png'));

fig6 = plotOrbit3D(mol, 'ShowPerigee', true);
saveas(fig6, fullfile(OUT_DIR, 'molniya_orbit3d.png'));

%% ── 4. GEO ───────────────────────────────────────────────────────────────────
geo = earthOrbit('geo');
fprintf('Geostationary Orbit\n');
fprintf('  a        = %.1f km\n', geo.a);
fprintf('  alt      = %.0f km\n', geo.alt_peri);
fprintf('  period   = %.4f h  (one sidereal day)\n\n', geo.period/3600);

fig7 = plotGroundTrack(geo, 'NumOrbits', 1, 'J2', false, 'ShowNodes', false);
saveas(fig7, fullfile(OUT_DIR, 'geo_groundtrack.png'));

fig8 = plotOrbit3D(geo, 'ShowNodeLine', false, 'ShowPerigee', false);
saveas(fig8, fullfile(OUT_DIR, 'geo_orbit3d.png'));

%% ── 5. TLE parsing example ───────────────────────────────────────────────────
% Nominal ISS TLE (synthetic, for demonstration)
line1 = '1 25544U 98067A   24001.50000000  .00010000  00000+0  18000-3 0  9997';
line2 = '2 25544  51.6434 260.4943 0001234  97.4120 262.7185 15.50000000000001';

tle = earthOrbit('tle', line1, line2);
fprintf('TLE orbit\n');
fprintf('  a        = %.1f km\n', tle.a);
fprintf('  i        = %.4f deg\n', tle.i);
fprintf('  RAAN     = %.4f deg\n', tle.RAAN);
fprintf('  e        = %.7f\n', tle.e);
fprintf('  period   = %.2f min\n\n', tle.period/60);

fig9 = plotGroundTrack(tle, 'NumOrbits', 3, 'J2', true);
saveas(fig9, fullfile(OUT_DIR, 'tle_groundtrack.png'));

%% ── 6. ECI <-> COE round-trip check ─────────────────────────────────────────
fprintf('Round-trip check: COE -> ECI -> COE\n');
a_in = 7000; e_in = 0.05; i_in = 45.0; RAAN_in = 60.0; om_in = 30.0; nu_in = 90.0;
[r, v] = coe2eci(a_in, e_in, i_in, RAAN_in, om_in, nu_in);
coe    = eci2coe(r, v);
fprintf('  a:     in=%.4f  out=%.4f  err=%.2e km\n',   a_in,    coe.a,    abs(coe.a - a_in));
fprintf('  e:     in=%.6f  out=%.6f  err=%.2e\n',       e_in,    coe.e,    abs(coe.e - e_in));
fprintf('  i:     in=%.4f  out=%.4f  err=%.2e deg\n',   i_in,    coe.i,    abs(coe.i - i_in));
fprintf('  RAAN:  in=%.4f  out=%.4f  err=%.2e deg\n',   RAAN_in, coe.RAAN, abs(coe.RAAN - RAAN_in));
fprintf('  omega: in=%.4f  out=%.4f  err=%.2e deg\n',   om_in,   coe.omega,abs(coe.omega - om_in));
fprintf('  nu:    in=%.4f  out=%.4f  err=%.2e deg\n\n', nu_in,   coe.nu,   abs(coe.nu - nu_in));

fprintf('Output images saved to: %s\n', OUT_DIR);
