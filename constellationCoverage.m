function cov = constellationCoverage(sats, duration_hr, varargin)
%CONSTELLATIONCOVERAGE  Combined coverage analysis for an entire constellation.
%
%   cov = constellationCoverage(sats, duration_hr)
%   cov = constellationCoverage(sats, duration_hr, Name, Value, ...)
%
%   Inputs:
%     sats        - 1xT struct array of earthOrbit structs
%     duration_hr - analysis duration (hours)
%
%   Options (Name-Value pairs):
%     'MinElevation' - minimum elevation angle (deg), default: 10
%     'GridRes'      - grid resolution (deg), default: 5
%     'StepSize'     - time step (s), default: 60
%
%   Returns struct compatible with plotCoverage:
%     coverage_frac    - Nlat x Nlon matrix, fraction of time covered
%     revisit_mean_hr  - Nlat x Nlon, mean gap between passes (NaN if <2)
%     revisit_max_hr   - Nlat x Nlon, max gap (NaN if never covered)
%     n_passes         - Nlat x Nlon, number of coverage events
%     lat_vec          - latitude grid vector (deg)
%     lon_vec          - longitude grid vector (deg)
%     duration_hr      - analysis duration (hours)
%     orb              - first satellite's orbit struct
%     min_elevation_deg
%     n_sats

p = inputParser;
addParameter(p, 'MinElevation', 10,  @isnumeric);
addParameter(p, 'GridRes',       5,  @isnumeric);
addParameter(p, 'StepSize',     60,  @isnumeric);
parse(p, varargin{:});
opts = p.Results;

R_E    = 6378.1363;      % km
om_E   = 7.2921150e-5;   % rad/s

duration_s = duration_hr * 3600;
dt         = opts.StepSize;
n_sats     = numel(sats);

%% ── Build lat/lon grid (matches coverageAnalysis exactly) ───────────────────
lat_vec = -90 : opts.GridRes : 90;    % 1xM
lon_vec = -180 : opts.GridRes : 180;  % 1xN_lon
M       = numel(lat_vec);
N_lon   = numel(lon_vec);
N_grid  = M * N_lon;

[LON, LAT] = meshgrid(lon_vec, lat_vec);   % MxN_lon each
LAT_v = LAT(:)';   % 1xN_grid
LON_v = LON(:)';

% Fixed ECEF target vectors (surface points)
r_tgt_ecef = R_E * [cosd(LAT_v) .* cosd(LON_v);
                     cosd(LAT_v) .* sind(LON_v);
                     sind(LAT_v)];   % 3xN_grid

%% ── Use epoch from first satellite for GMST ─────────────────────────────────
epoch_jd   = sats(1).epoch_jd;
GMST_epoch = mod(280.46061837 + 360.98564736629 * (epoch_jd - 2451545.0), 360);
om_E_deg   = rad2deg(om_E);   % deg/s

sin_MinElev = sind(opts.MinElevation);

%% ── Propagate all satellites ─────────────────────────────────────────────────
% Build uniform time vector
t_vec = (0 : dt : duration_s)';
if t_vec(end) < duration_s
    t_vec(end+1) = duration_s;
end
N_t = numel(t_vec);

% Propagate each satellite and store ECI positions (3 x N_t each)
r_eci_all = zeros(3, N_t, n_sats);

mu_E    = 398600.4418;
J2      = 1.08262668e-3;

for s = 1:n_sats
    orb_s   = sats(s);
    a_s     = orb_s.a;
    e_s     = orb_s.e;
    i_deg_s = orb_s.i;
    RAAN0_s = orb_s.RAAN;
    om0_s   = orb_s.omega;
    M0_s    = orb_s.M0;
    n_s     = 2*pi / orb_s.period;

    % J2 secular drift (matches propagateOrbit 'j2' method)
    p_orb    = a_s * (1 - e_s^2);
    factor   = 1.5 * n_s * J2 * (R_E / p_orb)^2;
    RAAN_dot = -factor * cosd(i_deg_s);
    om_dot   =  factor * (2.5 * cosd(i_deg_s)^2 - 1.0);
    n_corr   =  n_s + factor * sqrt(1 - e_s^2) * (1.5*cosd(i_deg_s)^2 - 0.5);

    for k = 1:N_t
        t_k   = t_vec(k);
        RAAN  = RAAN0_s + rad2deg(RAAN_dot) * t_k;
        omega = om0_s   + rad2deg(om_dot)   * t_k;
        M_rad = mod(deg2rad(M0_s) + n_corr * t_k, 2*pi);
        E     = keplerSolve(M_rad, e_s);
        nu    = 2 * atan2(sqrt(1+e_s)*sin(E/2), sqrt(1-e_s)*cos(E/2));
        [rv, ~] = coe2eci(a_s, e_s, i_deg_s, RAAN, omega, rad2deg(nu));
        r_eci_all(:, k, s) = rv;
    end
end

%% ── Statistics arrays (1 x N_grid) ─────────────────────────────────────────
visible_count = zeros(1, N_grid);
n_passes_vec  = zeros(1, N_grid);
in_cov        = false(1, N_grid);
ever_seen     = false(1, N_grid);
gap_t0        = zeros(1, N_grid);
total_gap_s   = zeros(1, N_grid);
max_gap_s     = zeros(1, N_grid);

%% ── Main time loop (vectorised over grid and satellites) ────────────────────
for k = 1:N_t
    t_k    = t_vec(k);
    GMST_k = mod(GMST_epoch + om_E_deg * t_k, 360);

    % Rotate ECEF targets -> ECI for this time step
    cG = cosd(GMST_k);
    sG = sind(GMST_k);
    r_tgt_eci = [ cG * r_tgt_ecef(1,:) - sG * r_tgt_ecef(2,:);
                  sG * r_tgt_ecef(1,:) + cG * r_tgt_ecef(2,:);
                       r_tgt_ecef(3,:)];   % 3xN_grid

    % Accumulate visibility across all satellites: any satellite above min elevation?
    visible = false(1, N_grid);

    for s = 1:n_sats
        r_sat = r_eci_all(:, k, s);   % 3x1

        % Slant vector and elevation
        delta_r  = r_sat - r_tgt_eci;   % 3xN_grid (broadcast)
        sin_el   = sum(r_tgt_eci .* delta_r, 1) ./ (R_E * vecnorm(delta_r, 2, 1));
        visible  = visible | (sin_el > sin_MinElev);
    end

    % Pass/gap tracking (same logic as coverageAnalysis)
    new_enter = visible & ~in_cov;
    new_exit  = ~visible & in_cov;

    n_passes_vec(new_enter) = n_passes_vec(new_enter) + 1;

    if k > 1 && any(new_exit)
        gap_t0(new_exit) = t_k;
    end

    if k > 1 && any(new_enter & ever_seen)
        reenter  = new_enter & ever_seen;
        gap_len  = t_k - gap_t0(reenter);
        total_gap_s(reenter) = total_gap_s(reenter) + gap_len;
        max_gap_s(reenter)   = max(max_gap_s(reenter), gap_len);
    end

    ever_seen(visible) = true;
    visible_count(visible) = visible_count(visible) + 1;
    in_cov = visible;
end

% Close gaps open at end of simulation
t_end = t_vec(end);
still_gapped = ever_seen & ~in_cov;
if any(still_gapped)
    gap_len_final = t_end - gap_t0(still_gapped);
    total_gap_s(still_gapped) = total_gap_s(still_gapped) + gap_len_final;
    max_gap_s(still_gapped)   = max(max_gap_s(still_gapped), gap_len_final);
end

%% ── Compute output quantities ───────────────────────────────────────────────
coverage_frac_v = visible_count / N_t;

% Mean revisit: total duration / number of passes (meaningful if >= 2 passes)
revisit_mean_v = nan(1, N_grid);
has2 = n_passes_vec >= 2;
revisit_mean_v(has2) = duration_s ./ n_passes_vec(has2) / 3600;

% Max gap
revisit_max_v = nan(1, N_grid);
revisit_max_v(ever_seen) = max_gap_s(ever_seen) / 3600;

%% ── Reshape to MxN_lon matrices ─────────────────────────────────────────────
coverage_frac   = reshape(coverage_frac_v,  M, N_lon);
revisit_mean_hr = reshape(revisit_mean_v,   M, N_lon);
revisit_max_hr  = reshape(revisit_max_v,    M, N_lon);
n_passes_mat    = reshape(n_passes_vec,     M, N_lon);

%% ── Output struct (matches coverageAnalysis format for plotCoverage) ─────────
cov = struct( ...
    'lat_vec',          lat_vec,              ...
    'lon_vec',          lon_vec,              ...
    'coverage_frac',    coverage_frac,        ...
    'revisit_mean_hr',  revisit_mean_hr,      ...
    'revisit_max_hr',   revisit_max_hr,       ...
    'n_passes',         n_passes_mat,         ...
    'duration_hr',      duration_hr,          ...
    'orb',              sats(1),              ...
    'min_elevation_deg',opts.MinElevation,    ...
    'n_sats',           n_sats);
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
