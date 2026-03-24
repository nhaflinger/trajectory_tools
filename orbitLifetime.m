function [result, varargout] = orbitLifetime(orb, varargin)
%ORBITLIFETIME  Estimate orbital lifetime due to atmospheric drag.
%
%   result = orbitLifetime(orb)
%   result = orbitLifetime(orb, Name, Value, ...)
%   [result, fig] = orbitLifetime(orb, ..., 'Plot', true)
%
%   Uses an orbit-averaged drag approach suitable for lifetime estimates over
%   months to years. Not a high-fidelity propagator.
%
%   Inputs:
%     orb  - orbit struct from earthOrbit()
%
%   Options (Name-Value pairs):
%     'CdAm'        - ballistic coefficient Cd*A/m in m^2/kg (default: 0.01)
%     'MaxYears'    - stop after this many years if reentry not reached (default: 30)
%     'StepOrbits'  - orbital periods per integration step (default: 10)
%     'Plot'        - generate figure (default: false)
%
%   Output struct fields:
%     t_days       - time history (days)
%     a_km         - semi-major axis history (km)
%     e            - eccentricity history
%     alt_peri_km  - perigee altitude history (km)
%     alt_apo_km   - apogee altitude history (km)
%     lifetime_days  - orbital lifetime (days)
%     lifetime_years - orbital lifetime (years)
%     reentered    - true if reentry criterion met
%     CdAm         - ballistic coefficient used

p = inputParser;
addParameter(p, 'CdAm',       0.01,  @isnumeric);
addParameter(p, 'MaxYears',   30,    @isnumeric);
addParameter(p, 'StepOrbits', 10,    @isnumeric);
addParameter(p, 'Plot',       false, @islogical);
parse(p, varargin{:});
opts = p.Results;

mu_SI  = 3.986004418e14;   % m^3/s^2
R_E_km = 6378.1363;        % km
R_E_m  = R_E_km * 1e3;     % m

a0_km  = orb.a;
e0     = orb.e;
CdAm   = opts.CdAm;

max_steps = ceil(opts.MaxYears * 365.25 * 86400 / (opts.StepOrbits * orb.period)) + 100;

% Pre-allocate history arrays
t_hist       = zeros(max_steps, 1);
a_hist       = zeros(max_steps, 1);
e_hist       = zeros(max_steps, 1);
alt_peri_hist = zeros(max_steps, 1);
alt_apo_hist  = zeros(max_steps, 1);

a_km  = a0_km;
e_cur = e0;
t_s   = 0.0;

n_stored  = 0;
reentered = false;

% 24 equally-spaced true anomaly sample points for orbit-averaging
N_nu  = 24;
nu_vec = linspace(0, 2*pi, N_nu + 1);
nu_vec = nu_vec(1:end-1);   % exclude 2*pi (same as 0)

max_t_s = opts.MaxYears * 365.25 * 86400;

step = 0;
while t_s <= max_t_s
    step = step + 1;

    % Current orbital elements
    a_m   = a_km * 1e3;                  % m
    p_m   = a_m * (1 - e_cur^2);        % semi-latus rectum, m
    h_m   = sqrt(mu_SI * p_m);          % specific angular momentum, m^2/s
    n_rad = sqrt(mu_SI / a_m^3);        % mean motion, rad/s
    T_orb = 2*pi / n_rad;               % period, s

    % Store current state
    n_stored = n_stored + 1;
    if n_stored > max_steps
        % Expand arrays if needed
        t_hist        = [t_hist;        zeros(max_steps,1)];  %#ok<AGROW>
        a_hist        = [a_hist;        zeros(max_steps,1)];  %#ok<AGROW>
        e_hist        = [e_hist;        zeros(max_steps,1)];  %#ok<AGROW>
        alt_peri_hist = [alt_peri_hist; zeros(max_steps,1)];  %#ok<AGROW>
        alt_apo_hist  = [alt_apo_hist;  zeros(max_steps,1)];  %#ok<AGROW>
    end

    alt_peri_km = a_km * (1 - e_cur) - R_E_km;
    alt_apo_km  = a_km * (1 + e_cur) - R_E_km;

    t_hist(n_stored)        = t_s / 86400;
    a_hist(n_stored)        = a_km;
    e_hist(n_stored)        = e_cur;
    alt_peri_hist(n_stored) = alt_peri_km;
    alt_apo_hist(n_stored)  = alt_apo_km;

    % Check reentry criterion: perigee below 80 km
    if alt_peri_km < 80
        reentered = true;
        break;
    end

    % Orbit-averaged density: weight samples by dt ∝ r^2 dν / h
    rho_weighted_sum = 0.0;
    weight_sum       = 0.0;

    for k = 1:N_nu
        nu_k = nu_vec(k);
        r_m  = p_m / (1 + e_cur * cos(nu_k));   % radius in m
        h_km_k = r_m / 1e3 - R_E_km;
        if h_km_k < 0
            h_km_k = 0;
        end
        rho_k  = atmDensityExp(h_km_k);
        % Weight: dt = r^2 dν / h  (Kepler's second law)
        wt_k   = r_m^2 / h_m;
        rho_weighted_sum = rho_weighted_sum + rho_k * wt_k;
        weight_sum       = weight_sum       + wt_k;
    end

    rho_mean = rho_weighted_sum / weight_sum;   % kg/m^3

    % Orbit-averaged da/dt (m/s)
    da_dt_m_s = -rho_mean * CdAm * sqrt(mu_SI * a_m);   % m/s

    % de/dt: King-Hele simplified — eccentricity decays proportionally
    % de/dt ≈ da/dt * e / (2*a)
    de_dt = (da_dt_m_s / (2 * a_m)) * e_cur;   % 1/s

    % Integration step: StepOrbits periods
    dt_step = opts.StepOrbits * T_orb;   % s

    % Forward Euler update
    da_km = (da_dt_m_s * dt_step) / 1e3;   % km
    de    = de_dt * dt_step;

    a_km  = a_km  + da_km;
    e_cur = e_cur + de;

    % Clamp eccentricity to [0, 0.999]
    e_cur = max(0, min(e_cur, 0.999));

    % Clamp a to not go below Earth radius + 10 km
    if a_km < R_E_km + 10
        a_km = R_E_km + 10;
        reentered = true;
        % Store final state
        n_stored = n_stored + 1;
        t_s = t_s + dt_step;
        alt_peri_km = a_km * (1 - e_cur) - R_E_km;
        alt_apo_km  = a_km * (1 + e_cur) - R_E_km;
        t_hist(n_stored)        = t_s / 86400;
        a_hist(n_stored)        = a_km;
        e_hist(n_stored)        = e_cur;
        alt_peri_hist(n_stored) = alt_peri_km;
        alt_apo_hist(n_stored)  = alt_apo_km;
        break;
    end

    t_s = t_s + dt_step;
end

% Trim arrays
t_days      = t_hist(1:n_stored);
a_km_out    = a_hist(1:n_stored);
e_out       = e_hist(1:n_stored);
alt_peri    = alt_peri_hist(1:n_stored);
alt_apo     = alt_apo_hist(1:n_stored);

lifetime_days  = t_days(end);
lifetime_years = lifetime_days / 365.25;

result = struct( ...
    't_days',        t_days,        ...
    'a_km',          a_km_out,      ...
    'e',             e_out,         ...
    'alt_peri_km',   alt_peri,      ...
    'alt_apo_km',    alt_apo,       ...
    'lifetime_days', lifetime_days, ...
    'lifetime_years',lifetime_years,...
    'reentered',     reentered,     ...
    'CdAm',          CdAm);

%% ── Print summary ──────────────────────────────────────────────────────────
fprintf('=== orbitLifetime ===\n');
fprintf('  Initial orbit : a=%.2f km, e=%.4f\n', orb.a, orb.e);
fprintf('  Initial alt   : perigee=%.1f km, apogee=%.1f km\n', orb.alt_peri, orb.alt_apo);
fprintf('  Cd*A/m        : %.4f m^2/kg\n', CdAm);
if reentered
    fprintf('  Lifetime      : %.1f days (%.2f years)  [REENTRY]\n', lifetime_days, lifetime_years);
else
    fprintf('  Lifetime      : > %.1f days (%.1f years)  [did not reenter]\n', lifetime_days, lifetime_years);
end
fprintf('=====================\n');

%% ── Optional plot ──────────────────────────────────────────────────────────
if opts.Plot
    bgCol  = [0.10 0.12 0.16];
    txtCol = [0.88 0.88 0.88];
    clrPeri = [0.30 0.75 0.93];
    clrApo  = [0.95 0.60 0.20];
    clrFill = [0.30 0.75 0.93];

    fig = figure('Color', bgCol, 'Position', [100 100 900 650]);

    % Top panel: altitude vs time
    ax1 = subplot(2,1,1, 'Parent', fig);
    set(ax1, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
        'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
    hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

    % Use log scale when altitude range spans more than one decade
    % (e.g. Molniya: 500 km perigee vs 40,000 km apogee)
    alt_ratio = max(alt_apo) / max(min(alt_peri(alt_peri > 0)), 1);
    use_log   = alt_ratio > 10;
    if use_log
        set(ax1, 'YScale', 'log');
    end

    % Only fill between perigee/apogee for nearly circular orbits — for highly
    % eccentric orbits the fill swamps the entire axes leaving nothing readable
    if orb.e < 0.1
        t_fill = [t_days; flipud(t_days)];
        y_fill = [alt_peri; flipud(alt_apo)];
        fill(ax1, t_fill, y_fill, clrFill, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    plot(ax1, t_days, alt_peri, '-', 'Color', clrPeri, 'LineWidth', 2.0, ...
        'DisplayName', 'Perigee alt');
    plot(ax1, t_days, alt_apo,  '-', 'Color', clrApo,  'LineWidth', 2.0, ...
        'DisplayName', 'Apogee alt');

    if reentered
        plot(ax1, t_days(end), alt_peri(end), 'v', 'MarkerSize', 10, ...
            'MarkerFaceColor', [0.90 0.30 0.30], 'MarkerEdgeColor', 'none', ...
            'DisplayName', 'Reentry');
    end

    yline(ax1, 80, '--', 'Color', [0.90 0.30 0.30], 'LineWidth', 1.2, ...
        'Label', '80 km reentry', 'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');

    xlabel(ax1, 'Time (days)', 'Color', txtCol);
    ylabel(ax1, 'Altitude (km)', 'Color', txtCol);

    orb_type = upper(orb.type);
    if strcmp(orb.type, 'circular')
        init_str = sprintf('%s %.0f km', orb_type, orb.alt_peri);
    else
        init_str = sprintf('%s peri=%.0f/apo=%.0f km', orb_type, orb.alt_peri, orb.alt_apo);
    end
    title(ax1, sprintf('Orbit Lifetime: %s, Cd*A/m=%.4f m^2/kg', init_str, CdAm), ...
        'Color', txtCol);

    leg1 = legend(ax1, 'Location', 'northeast');
    set(leg1, 'TextColor', txtCol, 'Color', bgCol+0.05, 'EdgeColor', [0.4 0.4 0.4]);

    % Bottom panel: eccentricity vs time
    ax2 = subplot(2,1,2, 'Parent', fig);
    set(ax2, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
        'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

    plot(ax2, t_days, e_out, '-', 'Color', [0.80 0.50 0.90], 'LineWidth', 2.0);

    if reentered
        plot(ax2, t_days(end), e_out(end), 'v', 'MarkerSize', 10, ...
            'MarkerFaceColor', [0.90 0.30 0.30], 'MarkerEdgeColor', 'none');
    end

    xlabel(ax2, 'Time (days)', 'Color', txtCol);
    ylabel(ax2, 'Eccentricity', 'Color', txtCol);
    title(ax2, 'Eccentricity vs Time', 'Color', txtCol);

    varargout{1} = fig;
else
    if nargout > 1
        varargout{1} = [];
    end
end

end

%% ── Local: USSA76 exponential atmosphere ────────────────────────────────────
function rho = atmDensityExp(h_km)
% ATMDENSITYEXP  USSA76 28-band exponential atmosphere model.
%   h_km - altitude above sea level (km)
%   rho  - atmospheric density (kg/m^3)

% Table: [alt_base(km), rho_base(kg/m^3), H(km)]
atm_table = [
      0,   1.225,       7.249;
     25,   3.899e-2,    6.349;
     30,   1.774e-2,    6.682;
     40,   3.972e-3,    7.554;
     50,   1.057e-3,    8.382;
     60,   3.206e-4,    7.714;
     70,   8.770e-5,    6.549;
     80,   1.905e-5,    5.799;
     90,   3.396e-6,    5.382;
    100,   5.297e-7,    5.877;
    110,   9.661e-8,    7.263;
    120,   2.438e-8,    9.473;
    130,   8.484e-9,   12.636;
    140,   3.845e-9,   16.149;
    150,   2.070e-9,   22.523;
    180,   5.464e-10,  29.740;
    200,   2.789e-10,  37.105;
    250,   7.248e-11,  45.546;
    300,   2.418e-11,  53.628;
    350,   9.158e-12,  53.298;
    400,   3.725e-12,  58.515;
    450,   1.585e-12,  60.828;
    500,   6.967e-13,  63.822;
    600,   1.454e-13,  71.835;
    700,   3.614e-14,  88.667;
    800,   1.170e-14, 124.640;
    900,   5.245e-15, 181.050;
   1000,   3.019e-15, 268.000];

h_km = max(h_km, 0);

n_bands = size(atm_table, 1);

% Find the correct band: use the last band whose base <= h_km
idx = find(atm_table(:,1) <= h_km, 1, 'last');
if isempty(idx)
    idx = 1;
end
if idx > n_bands
    idx = n_bands;
end

h_base = atm_table(idx, 1);
rho0   = atm_table(idx, 2);
H      = atm_table(idx, 3);

rho = rho0 * exp(-(h_km - h_base) / H);
end
