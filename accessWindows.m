function wins = accessWindows(orb, target_lat, target_lon, duration_s, varargin)
%ACCESSWINDOWS  Compute satellite visibility windows over a ground target.
%
%   wins = accessWindows(orb, target_lat, target_lon, duration_s)
%   wins = accessWindows(orb, target_lat, target_lon, duration_s, Name, Value, ...)
%
%   Inputs:
%     orb         - orbit struct returned by earthOrbit()
%     target_lat  - ground target geodetic latitude (deg)
%     target_lon  - ground target geodetic longitude (deg)
%     duration_s  - analysis duration (s)
%
%   Options (Name-Value pairs):
%     'MinElev'   - minimum elevation angle (deg), default 5
%     'Method'    - propagation method, default 'j2'
%     'StepSize'  - time step (s), default min(period/180, 30)
%
%   Returns a struct array. Each element has:
%     start_s      - window start time (s from epoch)
%     stop_s       - window stop time (s)
%     duration_s   - window duration (s)
%     max_elev_deg - maximum elevation angle during window (deg)
%     start_jd     - Julian Date of window start
%     stop_jd      - Julian Date of window stop

R_E  = 6378.1363;     % km
om_E = 7.2921150e-5;  % rad/s

p = inputParser;
addParameter(p, 'MinElev',  5,    @isnumeric);
addParameter(p, 'Method',  'j2', @ischar);
addParameter(p, 'StepSize', [],   @isnumeric);
parse(p, varargin{:});
opts = p.Results;

% Default step size: higher resolution than propagateOrbit default
if isempty(opts.StepSize)
    dt = min(orb.period / 180, 30);
else
    dt = opts.StepSize;
end

%% ── Propagate orbit ────────────────────────────────────────────────────────
traj = propagateOrbit(orb, duration_s, ...
    'Method',   opts.Method,   ...
    'StepSize', dt);

t_vec = traj.t;           % Nx1
r_eci = traj.r_eci;       % 3xN
N_t   = numel(t_vec);

%% ── Target ECEF position (fixed on Earth's surface) ───────────────────────
r_tgt_ecef = R_E * [cosd(target_lat)*cosd(target_lon);
                     cosd(target_lat)*sind(target_lon);
                     sind(target_lat)];

%% ── Compute GMST at epoch ──────────────────────────────────────────────────
GMST_epoch = mod(280.46061837 + 360.98564736629 * (orb.epoch_jd - 2451545.0), 360);
om_E_deg   = rad2deg(om_E);   % deg/s

%% ── Vectorised elevation computation ──────────────────────────────────────
% GMST at each time step (deg)
GMST_vec = mod(GMST_epoch + om_E_deg * t_vec, 360);   % Nx1

% Rotate ECEF target -> ECI at each time step
% x_eci =  cos(GMST)*x_ecef - sin(GMST)*y_ecef
% y_eci =  sin(GMST)*x_ecef + cos(GMST)*y_ecef
% z_eci =  z_ecef
cG = cosd(GMST_vec);   % Nx1
sG = sind(GMST_vec);   % Nx1

% r_tgt_eci: 3xN
r_tgt_eci = [ cG' * r_tgt_ecef(1) - sG' * r_tgt_ecef(2);
              sG' * r_tgt_ecef(1) + cG' * r_tgt_ecef(2);
              repmat(r_tgt_ecef(3), 1, N_t)];

% unit vector from target center outward (up direction)
up_hat = r_tgt_eci ./ vecnorm(r_tgt_eci, 2, 1);   % 3xN (constant magnitude = R_E)

% slant vector from target to satellite
delta_r = r_eci - r_tgt_eci;   % 3xN

% elevation angle
sin_el = sum(up_hat .* delta_r, 1) ./ vecnorm(delta_r, 2, 1);   % 1xN
elev   = asind(sin_el);                                           % 1xN (deg)

%% ── Detect windows (rising/falling edge) ──────────────────────────────────
above = elev > opts.MinElev;    % 1xN logical

% Rising edges: transition false->true; Falling: true->false
rising  = find(~above(1:end-1) &  above(2:end)) + 1;  % index where coverage begins
falling = find( above(1:end-1) & ~above(2:end));       % index where coverage ends

% Handle case where propagation starts/ends inside a window
if above(1)
    rising = [1, rising];
end
if above(end)
    falling = [falling, N_t];
end

% Pair up rising/falling edges
n_wins = min(numel(rising), numel(falling));

%% ── Build output struct ────────────────────────────────────────────────────
empty_win = struct('start_s', {}, 'stop_s', {}, 'duration_s', {}, ...
                   'max_elev_deg', {}, 'start_jd', {}, 'stop_jd', {});

if n_wins == 0
    wins = empty_win;
    return;
end

wins = empty_win;
for w = 1:n_wins
    i0 = rising(w);
    i1 = falling(w);
    if i0 > i1
        continue;
    end
    t0 = t_vec(i0);
    t1 = t_vec(i1);
    max_el = max(elev(i0:i1));

    wins(end+1).start_s      = t0;
    wins(end).stop_s         = t1;
    wins(end).duration_s     = t1 - t0;
    wins(end).max_elev_deg   = max_el;
    wins(end).start_jd       = orb.epoch_jd + t0 / 86400;
    wins(end).stop_jd        = orb.epoch_jd + t1 / 86400;
end
end
