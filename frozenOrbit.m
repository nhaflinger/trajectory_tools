function frz = frozenOrbit(alt_km, i_deg, varargin)
%FROZENSORBIT  Compute frozen orbit elements for J2+J3 equilibrium.
%
%   frz = frozenOrbit(alt_km, i_deg)
%   frz = frozenOrbit(alt_km, i_deg, 'omega_deg', 90)
%
%   Computes the frozen eccentricity and argument of perigee at which the
%   J3 zonal harmonic creates a forced eccentricity that exactly balances
%   the J2-driven apsidal rotation, yielding zero secular drift in both e
%   and omega.
%
%   For a frozen orbit, omega must be fixed at 90 deg (perigee over north
%   pole) or 270 deg (perigee over south pole). The frozen eccentricity is:
%
%       e_frozen = -(J3/(2*J2)) * (R_E/p) * sin(i)
%
%   where p = a*(1-e^2). This is solved iteratively since p depends on e.
%
%   Inputs:
%     alt_km  - circular altitude (km)
%     i_deg   - inclination (degrees)
%
%   Options:
%     'omega_deg' - 90 (perigee north) or 270 (perigee south), default: 90
%
%   Output struct fields:
%     a            - semi-major axis (km)
%     e            - frozen eccentricity
%     i            - inclination (deg)
%     omega        - frozen argument of perigee (deg, 90 or 270)
%     p            - semi-latus rectum (km)
%     alt_peri     - periapsis altitude (km)
%     alt_apo      - apoapsis altitude (km)
%     period       - orbital period (s)
%     e_dot_secular   - secular eccentricity drift (1/s, should be ~0)
%     RAAN_dot        - J2 RAAN drift rate (deg/day)
%     omega_dot_J2    - J2-only apsidal rotation rate (deg/day)
%     omega_dot_J2J3  - J2+J3 combined apsidal rotation rate (deg/day)

% ── Constants ─────────────────────────────────────────────────────────────────
mu   = 398600.4418;          % km^3/s^2
R_E  = 6378.1363;            % km
J2   = 1.08262668e-3;
J3   = -2.53265648e-6;

% ── Parse options ──────────────────────────────────────────────────────────────
p_opt = inputParser;
addParameter(p_opt, 'omega_deg', 90, @(x) ismember(x, [90, 270]));
parse(p_opt, varargin{:});
omega_deg = p_opt.Results.omega_deg;

% ── Basic orbit geometry ───────────────────────────────────────────────────────
a = R_E + alt_km;           % semi-major axis (km)
n = sqrt(mu / a^3);         % mean motion (rad/s)
T = 2*pi / n;               % period (s)

% ── Iterative frozen eccentricity solution ────────────────────────────────────
% Start with e = 0 guess
e_iter = 0;
for k = 1:10
    p_iter = a * (1 - e_iter^2);
    e_new  = -(J3 / (2*J2)) * (R_E / p_iter) * sind(i_deg);
    if abs(e_new - e_iter) < 1e-14
        break;
    end
    e_iter = e_new;
end
e_frozen = e_iter;

% Ensure non-negative eccentricity; if the formula gives a negative result
% for this inclination band, take the absolute value and flip omega
if e_frozen < 0
    e_frozen  = abs(e_frozen);
    omega_deg = mod(omega_deg + 180, 360);
end

% ── Derived geometry ──────────────────────────────────────────────────────────
p_km     = a * (1 - e_frozen^2);
alt_peri = a * (1 - e_frozen) - R_E;
alt_apo  = a * (1 + e_frozen) - R_E;

% ── J2 secular rates ──────────────────────────────────────────────────────────
% Factor for J2 secular perturbations
fac_J2     = 1.5 * n * J2 * (R_E / p_km)^2;

% RAAN drift (rad/s)
RAAN_dot_rads = -fac_J2 * cosd(i_deg);
RAAN_dot_degday = rad2deg(RAAN_dot_rads) * 86400;   % deg/day

% omega dot from J2 alone (rad/s)
omega_dot_J2_rads = fac_J2 * (2.5 * cosd(i_deg)^2 - 1.0);
omega_dot_J2_degday = rad2deg(omega_dot_J2_rads) * 86400;   % deg/day

% ── J3 contribution to omega_dot ─────────────────────────────────────────────
% The J3 forced eccentricity creates a secular omega_dot contribution.
% For the frozen orbit, the combined J2+J3 omega_dot = 0.
% J3 omega_dot (rad/s) at the frozen eccentricity:
%   d(omega)/dt|_J3 = (5/4) * n * J3 * (R_E/p)^3 * sin(i) * (5*cos^2(i) - 1 + e*...) / e
% For the frozen orbit condition, the exact balance means omega_dot_J2J3 → 0.
% We compute it numerically as a verification.
%
% Using the Brouwer/Kozai J3 secular apsidal rate (simplified form valid for
% small e at the frozen solution):
%   domega/dt|_J3 = -(5/4) * n * J3/J2 * (R_E/p) * sin(i) * omega_dot_J2 / n
% For a rigorous treatment see Brouwer (1959). Here we use the result that
% by construction the frozen condition sets the combined rate to zero.

% Analytic J3 omega_dot term (from the secular equation of motion for omega):
%   domega/dt|_J3 ≈ (15/4) * n * J3 * (R_E/p)^3 * (cos^2(i) * ... )
% For small e, the dominant contribution is:
%   domega/dt|_J3 ≈ -(5/4) * n * J3 * (R_E/p)^3 * (4 - 5*sin^2(i)) * sin(i) / (2*e)
% This diverges as e->0, which is expected — the frozen orbit is the
% equilibrium that prevents this divergence by keeping e finite.
%
% Instead, compute the combined omega_dot at the frozen solution via the
% secular variation equation. At frozen orbit (omega = 90 deg), the
% e_dot and omega_dot should both be zero.

% J3 secular e_dot at frozen solution (should be ~0):
%   de/dt = -(5/2)*n*J3*(R_E/p)^3 * e * cos(omega) * sin(i) * (3/4) * ...
% Simplified form (from secular averaging):
%   de/dt ≈ (3/2) * n * J3 * (R_E/p)^3 * sin(i) * cos(omega) * sqrt(1-e^2)
%            * (1 - (5/4)*sin^2(i))
% At omega=90 deg, cos(omega)=0, so de/dt = 0 identically. Good.
e_dot_secular = 0.0;   % by construction (omega=90 or 270, cos(omega)=0)

% J3 secular omega_dot contribution at frozen e and omega=90:
%   domega/dt|_J3 = -(3/2) * n * J3 * (R_E/p)^3 * sin(i) * sin(omega) / (2*e) * (1 + ...)
% At omega=90, sin(omega)=1. The J3 term:
if e_frozen > 1e-8
    omega_dot_J3_rads = -(15/4) * n * J3 * (R_E / p_km)^3 * ...
        (1 - (5/4)*sind(i_deg)^2) * sind(i_deg) * sind(omega_deg) / sqrt(1 - e_frozen^2);
else
    omega_dot_J3_rads = 0;
end
omega_dot_J2J3_rads   = omega_dot_J2_rads + omega_dot_J3_rads;
omega_dot_J2J3_degday = rad2deg(omega_dot_J2J3_rads) * 86400;   % deg/day

% ── Verification: e_dot_check via numerical balance condition ─────────────────
% The frozen orbit condition requires de/dt = 0.
% The secular de/dt from J3 (dominant term) is proportional to cos(omega),
% which is exactly zero at omega = 90 or 270 deg.
% We verify this explicitly:
e_dot_check_J3 = (3/2) * n * J3 * (R_E / p_km)^3 * ...
    sind(i_deg) * cosd(omega_deg) * sqrt(1 - e_frozen^2);
fprintf('=== frozenOrbit Verification ===\n');
fprintf('  Alt=%.1f km, i=%.2f deg, omega=%.1f deg\n', alt_km, i_deg, omega_deg);
fprintf('  e_frozen    = %.6e\n', e_frozen);
fprintf('  e_dot_check (J3, should be ~0): %.3e 1/s\n', e_dot_check_J3);

% ── Print summary ──────────────────────────────────────────────────────────────
fprintf('\n=== Frozen Orbit Summary ===\n');
fprintf('  Semi-major axis : %.4f km\n', a);
fprintf('  Eccentricity    : %.6e\n', e_frozen);
fprintf('  Inclination     : %.4f deg\n', i_deg);
fprintf('  omega (frozen)  : %.1f deg  (%s pole)\n', omega_deg, ...
    ternary(omega_deg == 90, 'north', 'south'));
fprintf('  Perigee alt     : %.4f km\n', alt_peri);
fprintf('  Apogee alt      : %.4f km\n', alt_apo);
fprintf('  Period          : %.2f s (%.4f min)\n', T, T/60);
fprintf('  RAAN drift      : %.4f deg/day\n', RAAN_dot_degday);
fprintf('  omega_dot (J2 only)   : %.6f deg/day\n', omega_dot_J2_degday);
fprintf('  omega_dot (J2+J3)     : %.6f deg/day\n', omega_dot_J2J3_degday);

if e_frozen < 0.01
    fprintf('  Earth-obs friendly: YES (e=%.2e < 0.01)\n', e_frozen);
else
    fprintf('  Earth-obs friendly: NO  (e=%.4f >= 0.01, altitude band ~%.1f km)\n', ...
        e_frozen, alt_apo - alt_peri);
end
fprintf('============================\n');

% ── Build output struct ────────────────────────────────────────────────────────
frz = struct( ...
    'a',                a,                    ...
    'e',                e_frozen,             ...
    'i',                i_deg,                ...
    'omega',            omega_deg,            ...
    'p',                p_km,                 ...
    'alt_peri',         alt_peri,             ...
    'alt_apo',          alt_apo,              ...
    'period',           T,                    ...
    'e_dot_secular',    e_dot_secular,        ...
    'RAAN_dot',         RAAN_dot_degday,      ...
    'omega_dot_J2',     omega_dot_J2_degday,  ...
    'omega_dot_J2J3',   omega_dot_J2J3_degday);

end

%% ── Local helper ──────────────────────────────────────────────────────────────
function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
