function result = betaAngle(orb, duration_days, varargin)
%BETAANGLE  Compute beta angle history and eclipse statistics for an Earth orbit.
%
%   result = betaAngle(orb, duration_days)
%   result = betaAngle(orb, duration_days, Name, Value, ...)
%
%   The beta angle is the angle between the orbital plane and the Sun direction:
%     beta = arcsin(dot(h_hat, sun_hat))
%
%   The RAAN precesses over time due to J2 secular drift. The Sun direction
%   also changes as Earth orbits the Sun.
%
%   Inputs:
%     orb           - orbit struct from earthOrbit()
%     duration_days - analysis duration (days)
%
%   Options:
%     'StepDays'    - time step in days (default: 1.0)
%
%   Output struct fields:
%     t_days         - time vector (days from epoch), Nx1
%     beta_deg       - beta angle (deg), Nx1
%     f_eclipse      - eclipse fraction per orbit (0-1), Nx1
%     t_eclipse_min  - eclipse duration per orbit (min), Nx1
%     RAAN_deg       - RAAN at each time step (deg), Nx1
%     beta_max_deg   - maximum beta angle (deg)
%     beta_min_deg   - minimum beta angle (deg)
%     eclipse_free   - true if orbit is always eclipse-free
%     orb            - copy of input orbit struct

%% ── Constants ────────────────────────────────────────────────────────────────
mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km
J2   = 1.08262668e-3;

%% ── Parse options ────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'StepDays', 1.0, @isnumeric);
parse(p, varargin{:});
step_days = p.Results.StepDays;

%% ── Time vector ──────────────────────────────────────────────────────────────
t_days = (0 : step_days : duration_days)';
if t_days(end) < duration_days
    t_days(end+1) = duration_days;
end
N = numel(t_days);

%% ── Orbit parameters ─────────────────────────────────────────────────────────
a    = orb.a;
e    = orb.e;
i    = orb.i;        % inclination (deg)
RAAN0 = orb.RAAN;   % initial RAAN (deg)

p_orb = a * (1 - e^2);
n     = sqrt(mu_E / a^3);   % rad/s

% J2 secular RAAN drift rate (rad/s)
RAAN_dot_rads = -1.5 * n * J2 * (R_E / p_orb)^2 * cosd(i);   % rad/s
RAAN_dot_degday = rad2deg(RAAN_dot_rads) * 86400;              % deg/day

% Earth's angular radius as seen from the orbit (for eclipse calculation)
rho_deg = asind(R_E / a);   % deg

%% ── Allocate outputs ─────────────────────────────────────────────────────────
beta_deg      = zeros(N, 1);
f_eclipse     = zeros(N, 1);
t_eclipse_min = zeros(N, 1);
RAAN_deg      = zeros(N, 1);

%% ── Main loop ────────────────────────────────────────────────────────────────
for k = 1:N
    t_d = t_days(k);

    % RAAN at this time (deg), accounting for J2 precession
    RAAN_k = RAAN0 + RAAN_dot_degday * t_d;
    RAAN_deg(k) = mod(RAAN_k, 360);

    % Orbit normal (angular momentum unit vector) in ECI
    % h_hat = [-sin(i)*sin(RAAN); sin(i)*cos(RAAN); cos(i)]
    h_hat = [-sind(i)*sind(RAAN_k); ...
              sind(i)*cosd(RAAN_k); ...
              cosd(i)];

    % Sun position at this time
    jd_now = orb.epoch_jd + t_d;
    [sun_hat, ~] = sunPosition(jd_now);

    % Beta angle
    beta_k = asind(dot(h_hat, sun_hat));
    beta_deg(k) = beta_k;

    % Eclipse fraction (cylindrical shadow model, spherical Earth)
    if abs(beta_k) < rho_deg
        f_e = (1/pi) * acos(sqrt(1 - (R_E/a)^2) / cosd(beta_k));
    else
        f_e = 0;
    end
    f_eclipse(k) = f_e;
    t_eclipse_min(k) = f_e * orb.period / 60;
end

%% ── Summary statistics ───────────────────────────────────────────────────────
beta_max_deg = max(beta_deg);
beta_min_deg = min(beta_deg);
eclipse_free = all(abs(beta_deg) >= rho_deg);

%% ── Print summary ────────────────────────────────────────────────────────────
fprintf('\n=== Beta Angle Analysis ===\n');
fprintf('  Orbit   : %s, a=%.1f km, i=%.2f deg, e=%.4f\n', ...
        orb.type, orb.a, orb.i, orb.e);
fprintf('  Duration: %.1f days, step=%.2f days\n', duration_days, step_days);
fprintf('  Beta max: %+.2f deg\n', beta_max_deg);
fprintf('  Beta min: %+.2f deg\n', beta_min_deg);
fprintf('  Eclipse boundary (rho): %.2f deg\n', rho_deg);
fprintf('  Max eclipse duration  : %.2f min/orbit\n', max(t_eclipse_min));
if eclipse_free
    fprintf('  Eclipse-free: YES (|beta| > rho at all times)\n');
else
    fprintf('  Eclipse-free: NO\n');
end

%% ── Build output struct ──────────────────────────────────────────────────────
result = struct( ...
    't_days',        t_days,        ...
    'beta_deg',      beta_deg,      ...
    'f_eclipse',     f_eclipse,     ...
    't_eclipse_min', t_eclipse_min, ...
    'RAAN_deg',      RAAN_deg,      ...
    'beta_max_deg',  beta_max_deg,  ...
    'beta_min_deg',  beta_min_deg,  ...
    'eclipse_free',  eclipse_free,  ...
    'rho_deg',       rho_deg,       ...
    'orb',           orb);
end
