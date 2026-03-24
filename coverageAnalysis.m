function cov = coverageAnalysis(orb, duration_s, varargin)
%COVERAGEANALYSIS  Compute coverage fraction and revisit time over a lat/lon grid.
%
%   cov = coverageAnalysis(orb, duration_s)
%   cov = coverageAnalysis(orb, duration_s, Name, Value, ...)
%
%   Inputs:
%     orb        - orbit struct returned by earthOrbit()
%     duration_s - analysis duration (s)
%
%   Options (Name-Value pairs):
%     'GridRes'   - grid resolution in degrees (default: 5)
%     'LatLim'    - [lat_min lat_max] (default: [-90 90])
%     'LonLim'    - [lon_min lon_max] (default: [-180 180])
%     'MinElev'   - minimum elevation angle (deg), default 5
%     'Method'    - propagation method, default 'j2'
%     'StepSize'  - time step (s), default min(period/180, 30)
%
%   Returns struct with:
%     lat_vec       - latitude grid vector (deg), 1xM
%     lon_vec       - longitude grid vector (deg), 1xN
%     coverage_frac - fraction of time with coverage, MxN
%     revisit_mean_hr - mean revisit time (hours), MxN (NaN where never covered)
%     revisit_max_hr  - maximum revisit gap (hours), MxN (NaN where never covered)
%     n_passes      - number of passes, MxN
%     duration_hr   - analysis duration in hours
%     orb           - input orbit struct

R_E  = 6378.1363;     % km
om_E = 7.2921150e-5;  % rad/s

p = inputParser;
addParameter(p, 'GridRes',  5,         @isnumeric);
addParameter(p, 'LatLim',   [-90 90],  @isnumeric);
addParameter(p, 'LonLim',   [-180 180],@isnumeric);
addParameter(p, 'MinElev',  5,         @isnumeric);
addParameter(p, 'Method',  'j2',       @ischar);
addParameter(p, 'StepSize', [],        @isnumeric);
parse(p, varargin{:});
opts = p.Results;

if isempty(opts.StepSize)
    dt = min(orb.period / 180, 30);
else
    dt = opts.StepSize;
end

%% ── Propagate orbit ────────────────────────────────────────────────────────
traj = propagateOrbit(orb, duration_s, ...
    'Method',   opts.Method, ...
    'StepSize', dt);

t_vec = traj.t;       % N_t x 1
r_eci = traj.r_eci;   % 3 x N_t
N_t   = numel(t_vec);

%% ── Build lat/lon grid ─────────────────────────────────────────────────────
lat_vec = opts.LatLim(1)  : opts.GridRes : opts.LatLim(2);   % 1xM
lon_vec = opts.LonLim(1)  : opts.GridRes : opts.LonLim(2);   % 1xN_lon
M       = numel(lat_vec);
N_lon   = numel(lon_vec);
N_grid  = M * N_lon;

% All grid points as ECEF vectors (3 x N_grid)
[LON, LAT] = meshgrid(lon_vec, lat_vec);   % MxN_lon each
LAT_v = LAT(:)';   % 1xN_grid
LON_v = LON(:)';

r_tgt_ecef = R_E * [cosd(LAT_v) .* cosd(LON_v);
                     cosd(LAT_v) .* sind(LON_v);
                     sind(LAT_v)];   % 3xN_grid

%% ── GMST setup ─────────────────────────────────────────────────────────────
GMST_epoch = mod(280.46061837 + 360.98564736629 * (orb.epoch_jd - 2451545.0), 360);
om_E_deg   = rad2deg(om_E);   % deg/s

sin_MinElev = sind(opts.MinElev);

%% ── Statistics arrays (1 x N_grid) ────────────────────────────────────────
visible_count = zeros(1, N_grid);   % total steps visible
n_passes_vec  = zeros(1, N_grid);   % number of coverage entries
in_cov        = false(1, N_grid);   % currently covered?
ever_seen     = false(1, N_grid);   % ever visible?
gap_t0        = zeros(1, N_grid);   % time when current gap started
total_gap_s   = zeros(1, N_grid);   % cumulative gap time
max_gap_s     = zeros(1, N_grid);   % max single gap duration

% Initialize gap start at t=0 for all points
gap_t0(:) = 0;

%% ── Main time loop (vectorised over grid) ──────────────────────────────────
for k = 1:N_t
    t_k    = t_vec(k);
    GMST_k = mod(GMST_epoch + om_E_deg * t_k, 360);

    % Rotate ECEF -> ECI for all grid points simultaneously
    cG = cosd(GMST_k);
    sG = sind(GMST_k);
    % Rz(-GMST): x_eci = cos(GMST)*x_ecef - sin(GMST)*y_ecef
    r_tgt_eci = [ cG * r_tgt_ecef(1,:) - sG * r_tgt_ecef(2,:);
                  sG * r_tgt_ecef(1,:) + cG * r_tgt_ecef(2,:);
                       r_tgt_ecef(3,:)];   % 3xN_grid

    % Slant vector: satellite - target
    delta_r = r_eci(:, k) - r_tgt_eci;   % 3xN_grid (broadcast)

    % sin(elevation) = dot(up_hat, delta_r) / norm(delta_r)
    % up_hat = r_tgt_eci / R_E  (all targets are on surface)
    sin_el  = sum(r_tgt_eci .* delta_r, 1) ./ (R_E * vecnorm(delta_r, 2, 1));   % 1xN_grid
    visible = sin_el > sin_MinElev;   % 1xN_grid logical

    % Points that just entered coverage (false -> true)
    new_enter = visible & ~in_cov;
    % Points that just left coverage (true -> false)
    new_exit  = ~visible & in_cov;

    % Update coverage entry count
    n_passes_vec(new_enter) = n_passes_vec(new_enter) + 1;

    % For points exiting coverage for the first time: no gap tracking
    % For points exiting coverage after first coverage:
    if k > 1 && any(new_exit)
        gap_t0(new_exit) = t_k;
    end

    % For points re-entering coverage: accumulate gap time
    if k > 1 && any(new_enter & ever_seen)
        reenter = new_enter & ever_seen;
        gap_len = t_k - gap_t0(reenter);
        total_gap_s(reenter) = total_gap_s(reenter) + gap_len;
        max_gap_s(reenter)   = max(max_gap_s(reenter), gap_len);
    end

    % First time seeing a point: record gap start from t=0 handled at init
    % Mark points as ever seen
    ever_seen(visible) = true;

    % Update visibility count and state
    visible_count(visible) = visible_count(visible) + 1;
    in_cov = visible;
end

% Close any gaps still open at end of simulation (for points currently not in coverage)
t_end = t_vec(end);
still_gapped = ever_seen & ~in_cov;
if any(still_gapped)
    gap_len_final = t_end - gap_t0(still_gapped);
    total_gap_s(still_gapped) = total_gap_s(still_gapped) + gap_len_final;
    max_gap_s(still_gapped)   = max(max_gap_s(still_gapped), gap_len_final);
end

%% ── Compute output quantities ──────────────────────────────────────────────
coverage_frac_v = visible_count / N_t;   % 1xN_grid

% Mean revisit: meaningful only if n_passes >= 2
revisit_mean_v = nan(1, N_grid);
has2 = n_passes_vec >= 2;
revisit_mean_v(has2) = duration_s ./ n_passes_vec(has2) / 3600;

% Max gap
revisit_max_v = nan(1, N_grid);
revisit_max_v(ever_seen) = max_gap_s(ever_seen) / 3600;

n_passes_grid_v = n_passes_vec;

%% ── Reshape to MxN_lon matrices ───────────────────────────────────────────
coverage_frac   = reshape(coverage_frac_v,   M, N_lon);
revisit_mean_hr = reshape(revisit_mean_v,     M, N_lon);
revisit_max_hr  = reshape(revisit_max_v,      M, N_lon);
n_passes_mat    = reshape(n_passes_grid_v,    M, N_lon);

%% ── Output struct ──────────────────────────────────────────────────────────
cov = struct( ...
    'lat_vec',        lat_vec,          ...
    'lon_vec',        lon_vec,          ...
    'coverage_frac',  coverage_frac,    ...
    'revisit_mean_hr',revisit_mean_hr,  ...
    'revisit_max_hr', revisit_max_hr,   ...
    'n_passes',       n_passes_mat,     ...
    'duration_hr',    duration_s/3600,  ...
    'orb',            orb);
end
