function wins = launchWindow(launch_lat, launch_lon, target_inc, target_RAAN, start_jd, varargin)
%LAUNCHWINDOW  Compute launch windows to reach a target orbit plane from a ground site.
%
%   wins = launchWindow(launch_lat, launch_lon, target_inc, target_RAAN, start_jd)
%   wins = launchWindow(..., Name, Value, ...)
%
%   launch_lat     - launch site geodetic latitude (deg, + North)
%   launch_lon     - launch site longitude (deg, East positive)
%   target_inc     - target orbit inclination (deg)
%   target_RAAN    - target RAAN at start_jd (deg)
%   start_jd       - start of search window (Julian Date)
%
%   Name-Value options:
%     'NDays'             - number of days to search (default: 1)
%     'SearchStep_s'      - time resolution in seconds (default: 60)
%     'RAANDrift_degPerDay' - secular RAAN drift of target orbit deg/day (default: 0)
%
%   Returns a struct array, one element per window found:
%     jd               - Julian Date of window
%     date_str         - formatted date/time string
%     type             - 'ascending' or 'descending'
%     azimuth_deg      - launch azimuth (deg from North, clockwise)
%     achieved_RAAN_deg
%     achieved_inc_deg
%     LMST_deg         - local mean sidereal time at launch

%% ── Parse options ────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'NDays',               1,   @isnumeric);
addParameter(p, 'SearchStep_s',        60,  @isnumeric);
addParameter(p, 'RAANDrift_degPerDay', 0,   @isnumeric);
parse(p, varargin{:});
n_days         = p.Results.NDays;
step_s         = p.Results.SearchStep_s;
raan_drift_dpd = p.Results.RAANDrift_degPerDay;

%% ── Validate inclination ─────────────────────────────────────────────────────
phi = launch_lat;   % geodetic latitude of launch site (deg)
if target_inc < abs(phi)
    error('launchWindow: target inclination (%.2f deg) is less than |launch latitude| (%.2f deg). Direct launch is not possible.', ...
          target_inc, abs(phi));
end

%% ── Launch azimuths ──────────────────────────────────────────────────────────
% Ascending node crossing launch azimuth
cos_az = cosd(target_inc) / cosd(phi);
cos_az = max(-1, min(1, cos_az));   % clamp for floating-point safety
az_asc  = asind(cos_az);            % 0..90 deg (prograde, heading NE)
az_desc = 180 - az_asc;             % heading SE

%% ── Argument of latitude of launch site on target orbit ─────────────────────
sin_uL = sind(phi) / sind(target_inc);
sin_uL = max(-1, min(1, sin_uL));
u_L = asind(sin_uL);                % 0..90 deg (ascending pass)

% Longitude difference from ascending node to sub-satellite point at u_L
delta_lon_asc  = atan2d(tand(u_L) * cosd(target_inc), 1);

% Descending pass: argument of latitude is in 2nd quadrant
u_L_desc       = 180 - u_L;
delta_lon_desc = atan2d(tand(u_L_desc) * cosd(target_inc), 1);

%% ── Time sweep ───────────────────────────────────────────────────────────────
step_days = step_s / 86400;
t_vec_jd  = (start_jd : step_days : start_jd + n_days)';
N         = numel(t_vec_jd);

% Pre-compute LMST at each time step
GMST_vec = mod(280.46061837 + 360.98564736629 * (t_vec_jd - 2451545.0), 360);
LMST_vec = mod(GMST_vec + launch_lon, 360);

% Pre-compute required LMST for ascending and descending windows
days_elapsed = t_vec_jd - start_jd;
RAAN_t       = target_RAAN + raan_drift_dpd * days_elapsed;   % target RAAN at each time

LMST_req_asc  = mod(RAAN_t + delta_lon_asc,  360);
LMST_req_desc = mod(RAAN_t + delta_lon_desc, 360);

%% ── Zero-crossing detection ──────────────────────────────────────────────────
% Wrap residuals into (-180, 180] to detect sign changes correctly
res_asc  = mod(LMST_vec - LMST_req_asc  + 180, 360) - 180;
res_desc = mod(LMST_vec - LMST_req_desc + 180, 360) - 180;

win_list = struct('jd', {}, 'date_str', {}, 'type', {}, ...
                  'azimuth_deg', {}, 'achieved_RAAN_deg', {}, ...
                  'achieved_inc_deg', {}, 'LMST_deg', {});

win_list = findCrossings(win_list, t_vec_jd, LMST_vec, RAAN_t, ...
                         res_asc,  az_asc,  target_inc, 'ascending');
win_list = findCrossings(win_list, t_vec_jd, LMST_vec, RAAN_t, ...
                         res_desc, az_desc, target_inc, 'descending');

% Sort by JD
if ~isempty(win_list)
    [~, sidx] = sort([win_list.jd]);
    win_list  = win_list(sidx);
end

wins = win_list;

%% ── Print table ──────────────────────────────────────────────────────────────
fprintf('\n=== Launch Windows ===\n');
fprintf('  Site  : lat=%.3f deg, lon=%.3f deg\n', launch_lat, launch_lon);
fprintf('  Target: i=%.2f deg, RAAN=%.2f deg\n', target_inc, target_RAAN);
fprintf('  Search: %.2f days from JD=%.4f\n', n_days, start_jd);
fprintf('  Found : %d window(s)\n\n', numel(wins));
if ~isempty(wins)
    fprintf('  %-4s  %-22s  %-11s  %-10s  %-12s\n', ...
            '#', 'Date (UTC)', 'Type', 'Az (deg)', 'RAAN (deg)');
    fprintf('  %s\n', repmat('-', 1, 68));
    for k = 1:numel(wins)
        w = wins(k);
        fprintf('  %-4d  %-22s  %-11s  %10.4f  %12.4f\n', ...
                k, w.date_str, w.type, w.azimuth_deg, w.achieved_RAAN_deg);
    end
end

end

%% ════════════════════════════════════════════════════════════════════════════
%%  Local helpers
%% ════════════════════════════════════════════════════════════════════════════

function win_list = findCrossings(win_list, t_jd, LMST_vec, RAAN_t, ...
                                  residual, azimuth, inc, pass_type)
% Find zero-crossings of residual (LMST - LMST_required) and refine by interp.
    for k = 1:numel(residual)-1
        r1 = residual(k);
        r2 = residual(k+1);
        % Zero crossing: sign change (ignore if both zero)
        if r1 * r2 < 0
            % Linear interpolation for precise crossing time
            frac = -r1 / (r2 - r1);
            jd_win  = t_jd(k)   + frac * (t_jd(k+1)   - t_jd(k));
            lmst_win = LMST_vec(k) + frac * (LMST_vec(k+1) - LMST_vec(k));
            raan_win = RAAN_t(k)   + frac * (RAAN_t(k+1)   - RAAN_t(k));

            [yr, mo, dy, hr, mn, sc] = jd_to_date(jd_win);
            date_str = sprintf('%04d-%02d-%02d %02d:%02d:%05.2f UTC', ...
                               yr, mo, dy, hr, mn, sc);

            entry.jd               = jd_win;
            entry.date_str         = date_str;
            entry.type             = pass_type;
            entry.azimuth_deg      = azimuth;
            entry.achieved_RAAN_deg = mod(raan_win, 360);
            entry.achieved_inc_deg  = inc;
            entry.LMST_deg          = mod(lmst_win, 360);

            win_list(end+1) = entry; %#ok<AGROW>
        end
    end
end

%% ── JD to calendar date ──────────────────────────────────────────────────────
function [yr, mo, dy, hr, mn, sc] = jd_to_date(jd)
    jd    = jd + 0.5;
    Z     = floor(jd);
    F     = jd - Z;
    if Z < 2299161
        A = Z;
    else
        alpha = floor((Z - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - floor(alpha/4);
    end
    B      = A + 1524;
    C      = floor((B - 122.1) / 365.25);
    D      = floor(365.25 * C);
    E      = floor((B - D) / 30.6001);
    dy_frac = B - D - floor(30.6001*E) + F;
    dy     = floor(dy_frac);
    hr_frac = (dy_frac - dy) * 24;
    hr     = floor(hr_frac);
    mn_frac = (hr_frac - hr) * 60;
    mn     = floor(mn_frac);
    sc     = (mn_frac - mn) * 60;
    if E < 14,  mo = E - 1;  else, mo = E - 13;  end
    if mo > 2,  yr = C - 4716;  else, yr = C - 4715;  end
end
