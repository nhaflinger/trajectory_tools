function result = propagateOrbit(orb, duration_s, varargin)
%PROPAGATEORBIT  Propagate an Earth orbit and return a time-history struct.
%
%   result = propagateOrbit(orb, duration_s)
%   result = propagateOrbit(orb, duration_s, Name, Value, ...)
%
%   Inputs:
%     orb        - orbit struct returned by earthOrbit()
%     duration_s - propagation duration in seconds
%
%   Options (Name-Value pairs):
%     'Method'     - 'kepler' (analytical 2-body), 'j2' (analytical J2 secular,
%                    default), 'numerical' (ode45 with J2 full periodic),
%                    'drag' (ode45 with J2 + exponential drag),
%                    'j3' (ode45 with J2+J3),
%                    'j4' (ode45 with J2+J3+J4),
%                    'srp' (ode45 with J2 + solar radiation pressure),
%                    'full' (ode45 with J2+J3+J4+drag+SRP)
%     'StepSize'   - output time step in seconds (default: min(period/360, 60))
%     'CdAm'       - Cd*A/m drag coefficient in m^2/kg (default: 0.01)
%     'CR_Am'      - reflectivity * area/mass in m^2/kg for SRP (default: 0.015)
%     'OutputCOE'  - also compute osculating COE at each step (default: false)
%
%   Output struct fields:
%     t       - time vector (s), Nx1
%     r_eci   - ECI position (km), 3xN
%     v_eci   - ECI velocity (km/s), 3xN
%     lat     - geodetic latitude (deg), Nx1
%     lon     - geodetic longitude (deg), Nx1
%     alt     - altitude above R_E (km), Nx1
%     orb     - copy of input orbit struct
%
%   If OutputCOE=true, also:
%     a, e, i_deg, RAAN, omega, M  (each Nx1)

p = inputParser;
addParameter(p, 'Method',    'j2',   @ischar);
addParameter(p, 'StepSize',  [],     @isnumeric);
addParameter(p, 'CdAm',      0.01,   @isnumeric);
addParameter(p, 'CR_Am',     0.015,  @isnumeric);
addParameter(p, 'OutputCOE', false,  @islogical);
parse(p, varargin{:});
opts = p.Results;

mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km
J2   = 1.08262668e-3;
J3   = -2.53265648e-6;
J4   = -1.61962159e-6;
om_E = 7.2921150e-5;  % rad/s

a     = orb.a;
e     = orb.e;
i_deg = orb.i;
RAAN0 = orb.RAAN;
om0   = orb.omega;
M0    = orb.M0;
T     = orb.period;
n     = 2*pi / T;

% Default step size
if isempty(opts.StepSize)
    dt = min(T / 360, 60);
else
    dt = opts.StepSize;
end

% Build time vector
t_vec = (0 : dt : duration_s)';
if t_vec(end) < duration_s
    t_vec(end+1) = duration_s;
end
N_t = numel(t_vec);

% GMST at epoch (deg)
GMST_epoch = mod(280.46061837 + 360.98564736629 * (orb.epoch_jd - 2451545.0), 360);
om_E_deg   = rad2deg(om_E);   % deg/s

method = lower(strtrim(opts.Method));
valid_methods = {'kepler','j2','numerical','drag','j3','j4','srp','full'};
if ~any(strcmp(method, valid_methods))
    error('propagateOrbit: unknown Method ''%s''. Valid: %s', method, strjoin(valid_methods, ', '));
end

%% ── Analytical methods (kepler / j2) ──────────────────────────────────────
if strcmp(method, 'kepler') || strcmp(method, 'j2')

    if strcmp(method, 'j2')
        p_orb    = a * (1 - e^2);
        factor   = 1.5 * n * J2 * (R_E / p_orb)^2;
        RAAN_dot = -factor * cosd(i_deg);                                % rad/s
        om_dot   =  factor * (2.5 * cosd(i_deg)^2 - 1.0);               % rad/s
        n_corr   =  n + factor * sqrt(1 - e^2) * (1.5*cosd(i_deg)^2 - 0.5);
    else
        RAAN_dot = 0;  om_dot = 0;  n_corr = n;
    end

    r_eci = zeros(3, N_t);
    v_eci = zeros(3, N_t);
    lat   = zeros(N_t, 1);
    lon   = zeros(N_t, 1);
    alt   = zeros(N_t, 1);

    if opts.OutputCOE
        a_out     = zeros(N_t, 1);
        e_out     = zeros(N_t, 1);
        i_out     = zeros(N_t, 1);
        RAAN_out  = zeros(N_t, 1);
        omega_out = zeros(N_t, 1);
        M_out     = zeros(N_t, 1);
    end

    for k = 1:N_t
        t     = t_vec(k);
        RAAN  = RAAN0 + rad2deg(RAAN_dot) * t;
        omega = om0   + rad2deg(om_dot)   * t;
        M_rad = mod(deg2rad(M0) + n_corr * t, 2*pi);
        E     = keplerSolve(M_rad, e);
        nu    = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));

        [rv, vv] = coe2eci(a, e, i_deg, RAAN, omega, rad2deg(nu));
        r_eci(:, k) = rv;
        v_eci(:, k) = vv;

        % ECI -> ECEF -> geodetic
        GMST_k  = mod(GMST_epoch + om_E_deg * t, 360);
        r_ecef  = eci2ecef(rv, GMST_k);
        lat(k)  = asind(r_ecef(3) / norm(r_ecef));
        lon(k)  = atan2d(r_ecef(2), r_ecef(1));
        alt(k)  = norm(r_ecef) - R_E;

        if opts.OutputCOE
            a_out(k)     = a;
            e_out(k)     = e;
            i_out(k)     = i_deg;
            RAAN_out(k)  = mod(RAAN, 360);
            omega_out(k) = mod(omega, 360);
            M_out(k)     = rad2deg(mod(M_rad, 2*pi));
        end
    end

%% ── Numerical methods ──────────────────────────────────────────────────────
elseif any(strcmp(method, {'numerical','drag','j3','j4','srp','full'}))

    % Initial ECI state from orbit struct
    y0 = [orb.r_vec(:); orb.v_vec(:)];

    % Build params struct for eom
    params.mu    = mu_E;
    params.R_E   = R_E;
    params.J2    = J2;
    params.J3    = J3;
    params.J4    = J4;
    params.om_E  = om_E;
    params.CdAm  = opts.CdAm;
    params.CR_Am = opts.CR_Am;
    params.jd_epoch = orb.epoch_jd;

    % Perturbation flags
    params.useDrag = strcmp(method, 'drag') || strcmp(method, 'full');
    params.useJ3   = any(strcmp(method, {'j3','j4','full'}));
    params.useJ4   = any(strcmp(method, {'j4','full'}));
    params.useSRP  = any(strcmp(method, {'srp','full'}));

    ode_opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    [t_out, y_out] = ode45(@(t,y) eomFull(t, y, params), ...
                           [0, duration_s], y0, ode_opts);

    % Interpolate to uniform t_vec
    % interp1(x, V, xq): x is N_steps x 1, V is N_steps x 3, xq is N_t x 1
    % result is N_t x 3; transpose to 3 x N_t
    r_eci_raw = interp1(t_out, y_out(:,1:3), t_vec, 'linear')';   % 3xN_t
    v_eci_raw = interp1(t_out, y_out(:,4:6), t_vec, 'linear')';   % 3xN_t

    r_eci = r_eci_raw;
    v_eci = v_eci_raw;

    lat   = zeros(N_t, 1);
    lon   = zeros(N_t, 1);
    alt   = zeros(N_t, 1);

    if opts.OutputCOE
        a_out     = zeros(N_t, 1);
        e_out     = zeros(N_t, 1);
        i_out     = zeros(N_t, 1);
        RAAN_out  = zeros(N_t, 1);
        omega_out = zeros(N_t, 1);
        M_out     = zeros(N_t, 1);
    end

    for k = 1:N_t
        GMST_k  = mod(GMST_epoch + om_E_deg * t_vec(k), 360);
        r_ecef  = eci2ecef(r_eci(:,k), GMST_k);
        lat(k)  = asind(r_ecef(3) / norm(r_ecef));
        lon(k)  = atan2d(r_ecef(2), r_ecef(1));
        alt(k)  = norm(r_ecef) - R_E;

        if opts.OutputCOE
            coe_k        = eci2coe(r_eci(:,k), v_eci(:,k));
            a_out(k)     = coe_k.a;
            e_out(k)     = coe_k.e;
            i_out(k)     = coe_k.i;
            RAAN_out(k)  = coe_k.RAAN;
            omega_out(k) = coe_k.omega;
            M_out(k)     = coe_k.M;
        end
    end

else
    % This branch should not be reached due to early validation above
    error('propagateOrbit: unknown Method ''%s''.', method);
end

%% ── Build output struct ────────────────────────────────────────────────────
result = struct( ...
    't',     t_vec,  ...
    'r_eci', r_eci,  ...
    'v_eci', v_eci,  ...
    'lat',   lat,    ...
    'lon',   lon,    ...
    'alt',   alt,    ...
    'orb',   orb);

if opts.OutputCOE
    result.a     = a_out;
    result.e     = e_out;
    result.i_deg = i_out;
    result.RAAN  = RAAN_out;
    result.omega = omega_out;
    result.M     = M_out;
end
end

%% ── Local: full equations of motion (J2+J3+J4+drag+SRP) ────────────────────
function dydt = eomFull(t, y, params)
    r = y(1:3);
    v = y(4:6);
    r_n = norm(r);

    mu  = params.mu;
    R_E = params.R_E;
    J2  = params.J2;

    % Two-body
    a = -mu / r_n^3 * r;

    % J2 perturbation (always on for numerical methods)
    z2   = (r(3) / r_n)^2;
    c2   = 1.5 * J2 * mu * R_E^2 / r_n^5;
    a_J2 = c2 * [r(1)*(5*z2 - 1); ...
                  r(2)*(5*z2 - 1); ...
                  r(3)*(5*z2 - 3)];
    a = a + a_J2;

    % J3 perturbation
    if params.useJ3
        s    = r(3) / r_n;
        c3   = (5/2) * mu * params.J3 * R_E^3 / r_n^6;  %#ok<PFBNS>
        a_J3 = c3 * [r(1)*s*(7*s^2 - 3); ...
                      r(2)*s*(7*s^2 - 3); ...
                      r_n * (3 - 30*s^2 + 35*s^4) / 5];
        a = a + a_J3;
    end

    % J4 perturbation
    if params.useJ4
        s    = r(3) / r_n;
        c4   = (15/8) * mu * params.J4 * R_E^4 / r_n^7;
        a_J4 = c4 * [r(1)*(21*s^4 - 14*s^2 + 1); ...
                      r(2)*(21*s^4 - 14*s^2 + 1); ...
                      r(3)*(63*s^4 - 70*s^2 + 15) / 3];
        a = a + a_J4;
    end

    % Drag perturbation
    if params.useDrag
        h_km  = r_n - R_E;
        rho   = atmDensityExp(h_km);
        v_rel = v - cross([0; 0; params.om_E], r);
        vr    = norm(v_rel);
        a_drag = -0.5 * rho * params.CdAm * vr * v_rel * 1e3;
        a = a + a_drag;
    end

    % Solar radiation pressure perturbation (cannonball model)
    if params.useSRP
        jd_now = params.jd_epoch + t / 86400;
        [sun_hat, ~] = sunPosition(jd_now);
        % P_srp [N/m^2] * CR_Am [m^2/kg] = [m/s^2] -> [km/s^2] = *1e-3
        P_srp  = 4.56e-6 * 1e-3;   % km/s^2 per (m^2/kg)
        a_srp  = -P_srp * params.CR_Am * sun_hat;

        % Cylindrical shadow check
        r_sun_dot = dot(r, sun_hat);
        r_perp    = norm(r - r_sun_dot * sun_hat);
        if r_sun_dot < 0 && r_perp < R_E
            a_srp = zeros(3, 1);   % in shadow
        end
        a = a + a_srp;
    end

    dydt = [v; a];
end

%% ── Local: equations of motion (legacy interface, kept for compatibility) ───
function dydt = eom(~, y, mu, R_E, J2, om_E, CdAm, do_drag)
    r = y(1:3);
    v = y(4:6);
    r_norm = norm(r);

    % Two-body
    a_grav = -mu / r_norm^3 * r;

    % J2 perturbation
    z2  = (r(3) / r_norm)^2;
    c   = 1.5 * J2 * mu * R_E^2 / r_norm^5;
    a_J2 = c * [r(1)*(5*z2 - 1);
                r(2)*(5*z2 - 1);
                r(3)*(5*z2 - 3)];

    a_total = a_grav + a_J2;

    % Drag perturbation
    if do_drag
        h_km  = r_norm - R_E;
        rho   = atmDensityExp(h_km);
        v_rel = v - cross([0; 0; om_E], r);   % km/s
        vr    = norm(v_rel);
        % rho [kg/m^3] * CdAm [m^2/kg] * v_rel^2 [km^2/s^2] -> need km/s^2
        % 0.5 * rho * CdAm * vr * v_rel  [kg/m^3 * m^2/kg * km/s * km/s]
        % = [m^-1 * km^2/s^2] = [1e3 km/s^2 / km] -> x 1e3 to get km/s^2
        a_drag = -0.5 * rho * CdAm * vr * v_rel * 1e3;
        a_total = a_total + a_drag;
    end

    dydt = [v; a_total];
end

%% ── Local: ECI -> ECEF rotation ─────────────────────────────────────────────
function r_ecef = eci2ecef(r_eci, GMST_deg)
    cG = cosd(GMST_deg);
    sG = sind(GMST_deg);
    r_ecef = [ cG*r_eci(1) + sG*r_eci(2);
              -sG*r_eci(1) + cG*r_eci(2);
               r_eci(3)];
end

%% ── Local: Kepler equation solver ───────────────────────────────────────────
function E = keplerSolve(M, e)
    E = M;
    for k = 1:50
        dE = (M - E + e*sin(E)) / (1 - e*cos(E));
        E  = E + dE;
        if abs(dE) < 1e-13, break; end
    end
end

%% ── Local: USSA76 exponential atmosphere ────────────────────────────────────
function rho = atmDensityExp(h_km)
    % Clamp to sea level minimum
    h_km = max(h_km, 0);

    % Table: [h_low, h_high, rho0 (kg/m^3), H (km)]
    tbl = [ ...
        0,    25,   1.225e+0,   7.249;
        25,   30,   3.899e-2,   6.349;
        30,   40,   1.774e-2,   6.682;
        40,   50,   3.972e-3,   7.554;
        50,   60,   1.057e-3,   8.382;
        60,   70,   3.206e-4,   7.714;
        70,   80,   8.770e-5,   6.549;
        80,   90,   1.905e-5,   5.799;
        90,  100,   3.396e-6,   5.382;
       100,  110,   5.297e-7,   5.877;
       110,  120,   9.661e-8,   7.263;
       120,  130,   2.438e-8,   9.473;
       130,  140,   8.484e-9,  12.636;
       140,  150,   3.845e-9,  16.149;
       150,  180,   2.070e-9,  22.523;
       180,  200,   5.464e-10, 29.740;
       200,  250,   2.789e-10, 37.105;
       250,  300,   7.248e-11, 45.546;
       300,  350,   2.418e-11, 53.628;
       350,  400,   9.518e-12, 53.298;
       400,  450,   3.725e-12, 58.515;
       450,  500,   1.585e-12, 60.828;
       500,  600,   6.967e-13, 63.822;
       600,  700,   1.454e-13, 71.835;
       700,  800,   3.614e-14, 88.667;
       800,  900,   1.170e-14, 124.64;
       900, 1000,   5.245e-15, 181.05;
      1000,  Inf,   3.019e-15, 268.00];

    % Find the correct band
    idx = find(h_km >= tbl(:,1) & h_km < tbl(:,2), 1);
    if isempty(idx)
        idx = size(tbl, 1);   % above 1000 km: use last row
    end
    rho0  = tbl(idx, 3);
    H     = tbl(idx, 4);
    h_low = tbl(idx, 1);
    rho   = rho0 * exp(-(h_km - h_low) / H);
end
