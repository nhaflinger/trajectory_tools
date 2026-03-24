function [result, fig] = dogLegTrade(launch_lat, launch_lon, target_inc, target_RAAN, start_jd, varargin)
%DOGLEGTRADE  Analyze the delta-V cost of launching at non-optimal times.
%
%   result = dogLegTrade(launch_lat, launch_lon, target_inc, target_RAAN, start_jd)
%   result = dogLegTrade(..., Name, Value, ...)
%   [result, fig] = dogLegTrade(..., 'Plot', true)
%
%   The function sweeps over a day at 1-minute resolution and computes the
%   dog-leg delta-V penalty for launching at each time compared to waiting
%   for the optimal (zero-cost) launch window.
%
%   Physics: At non-optimal launch times the achieved RAAN differs from the
%   target. A plane change after insertion corrects this at a delta-V cost:
%     dv_dogleg = 2 * v_insertion * sin(delta_RAAN / 2)
%   This is conservative — an in-flight dog leg during ascent is cheaper.
%
%   Inputs:
%     launch_lat   - launch site geodetic latitude (deg, +North)
%     launch_lon   - launch site longitude (deg, +East)
%     target_inc   - target orbit inclination (deg)
%     target_RAAN  - target RAAN (deg)
%     start_jd     - analysis start epoch (Julian Date)
%
%   Options:
%     'NDays'           - analysis period in days (default: 1)
%     'InsertionAlt_km' - orbit altitude for velocity (default: 300 km)
%     'Plot'            - generate figure (default: false)
%
%   Output struct fields:
%     t_hr            - time from start (hours), Nx1
%     jd              - Julian Date, Nx1
%     dv_dogleg_m_s   - dog-leg delta-V penalty (m/s), Nx1
%     RAAN_asc        - RAAN from ascending launch, Nx1
%     RAAN_desc       - RAAN from descending launch, Nx1
%     window_times_hr - times of zero-cost windows (hours)
%     min_dv_m_s      - minimum DV over period (m/s)
%     max_dv_m_s      - maximum DV over period (m/s)

%% ── Constants ────────────────────────────────────────────────────────────────
mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km

%% ── Parse options ────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'NDays',           1,     @isnumeric);
addParameter(p, 'InsertionAlt_km', 300,   @isnumeric);
addParameter(p, 'Plot',            false, @islogical);
parse(p, varargin{:});
n_days      = p.Results.NDays;
alt_ins     = p.Results.InsertionAlt_km;
do_plot     = p.Results.Plot;

%% ── Derived parameters ───────────────────────────────────────────────────────
a_ins   = R_E + alt_ins;
v_ins   = sqrt(mu_E / a_ins);   % km/s (circular velocity at insertion altitude)

phi = launch_lat;

% Validate inclination
if target_inc < abs(phi)
    error('dogLegTrade: target inclination (%.2f deg) less than |launch lat| (%.2f deg)', ...
          target_inc, abs(phi));
end

%% ── Launch geometry (same as launchWindow.m) ────────────────────────────────
% Ascending node crossing launch azimuth
cos_az = cosd(target_inc) / cosd(phi);
cos_az = max(-1, min(1, cos_az));

% Argument of latitude of launch site on target orbit
sin_uL = sind(phi) / sind(target_inc);
sin_uL = max(-1, min(1, sin_uL));
u_L    = asind(sin_uL);   % 0..90 deg (ascending pass)

% Longitude offset from ascending node to sub-satellite point
delta_lon_asc  = atan2d(tand(u_L) * cosd(target_inc), 1);

% Descending pass: argument of latitude is in 2nd quadrant
u_L_desc       = 180 - u_L;
delta_lon_desc = atan2d(tand(u_L_desc) * cosd(target_inc), 1);

%% ── Time sweep ───────────────────────────────────────────────────────────────
step_s    = 60;   % 1 minute resolution
step_days = step_s / 86400;
jd_vec    = (start_jd : step_days : start_jd + n_days)';
N         = numel(jd_vec);
t_hr      = (jd_vec - start_jd) * 24;   % hours from start

% GMST and LMST at each time step
GMST_vec = mod(280.46061837 + 360.98564736629 * (jd_vec - 2451545.0), 360);
LMST_vec = mod(GMST_vec + launch_lon, 360);

% RAAN achieved by ascending/descending launch at each time
RAAN_asc  = mod(LMST_vec - delta_lon_asc,  360);
RAAN_desc = mod(LMST_vec - delta_lon_desc, 360);

%% ── Dog-leg DV calculation ───────────────────────────────────────────────────
dv_dogleg_m_s = zeros(N, 1);

for k = 1:N
    % Angular difference between target RAAN and best achievable RAAN
    dRAAN_asc  = mod(RAAN_asc(k)  - target_RAAN + 180, 360) - 180;
    dRAAN_desc = mod(RAAN_desc(k) - target_RAAN + 180, 360) - 180;

    % Pick the option with smallest RAAN error
    if abs(dRAAN_asc) <= abs(dRAAN_desc)
        best_dRAAN = dRAAN_asc;
    else
        best_dRAAN = dRAAN_desc;
    end

    % Plane change DV (upper bound)
    dv_km_s = 2 * v_ins * abs(sind(best_dRAAN / 2));
    dv_dogleg_m_s(k) = dv_km_s * 1e3;   % m/s
end

%% ── Find natural launch windows (zero-cost) ─────────────────────────────────
thr_m_s = 1.0;   % threshold below which we consider it a "window" (m/s)
is_window = dv_dogleg_m_s < thr_m_s;
% Find transitions into window
win_starts = find(diff([0; is_window]) > 0);
window_times_hr = t_hr(win_starts);

%% ── Summary statistics ───────────────────────────────────────────────────────
min_dv_m_s = min(dv_dogleg_m_s);
max_dv_m_s = max(dv_dogleg_m_s);

%% ── Print summary ────────────────────────────────────────────────────────────
fprintf('\n=== Dog-Leg Trade Analysis ===\n');
fprintf('  Launch site   : lat=%.2f deg, lon=%.2f deg\n', launch_lat, launch_lon);
fprintf('  Target orbit  : i=%.1f deg, RAAN=%.1f deg\n', target_inc, target_RAAN);
fprintf('  Insertion alt : %.0f km (v=%.4f km/s)\n', alt_ins, v_ins);
fprintf('  Analysis period: %.1f day(s) from JD=%.4f\n', n_days, start_jd);
fprintf('  Max dog-leg DV: %.1f m/s\n', max_dv_m_s);
fprintf('  Min dog-leg DV: %.2f m/s\n', min_dv_m_s);
fprintf('  Natural windows (<%.0f m/s): %d found\n', thr_m_s, numel(window_times_hr));
for k = 1:min(numel(window_times_hr), 10)
    fprintf('    Window %d: t = %.2f hr from start\n', k, window_times_hr(k));
end

%% ── Build output struct ──────────────────────────────────────────────────────
result = struct( ...
    't_hr',            t_hr,            ...
    'jd',              jd_vec,          ...
    'dv_dogleg_m_s',   dv_dogleg_m_s,   ...
    'RAAN_asc',        RAAN_asc,        ...
    'RAAN_desc',       RAAN_desc,       ...
    'window_times_hr', window_times_hr, ...
    'min_dv_m_s',      min_dv_m_s,     ...
    'max_dv_m_s',      max_dv_m_s,     ...
    'target_RAAN',     target_RAAN,    ...
    'target_inc',      target_inc,     ...
    'launch_lat',      launch_lat,     ...
    'launch_lon',      launch_lon,     ...
    'alt_ins_km',      alt_ins,        ...
    'v_ins_km_s',      v_ins);

%% ── Optional plot ────────────────────────────────────────────────────────────
fig = [];
if do_plot
    bgCol  = [0.10 0.12 0.16];
    txtCol = [0.88 0.88 0.88];
    cyanC  = [0.30 0.75 0.93];
    orngC  = [0.95 0.60 0.20];
    greenC = [0.30 0.90 0.40];

    fig = figure('Color', bgCol, 'Position', [100 100 900 600]);

    title_str = sprintf('Dog-Leg Trade | site: %.1f°N %.1f°E | target: i=%.1f° RAAN=%.1f°', ...
        launch_lat, launch_lon, target_inc, target_RAAN);

    %% Top panel: DV penalty vs time
    ax1 = subplot(2, 1, 1);
    set(ax1, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
             'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
    hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

    % Fill area under DV curve
    patch(ax1, [t_hr; flipud(t_hr)], [dv_dogleg_m_s; zeros(size(dv_dogleg_m_s))], ...
          orngC, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(ax1, t_hr, dv_dogleg_m_s, '-', 'Color', orngC, 'LineWidth', 1.5);

    % Mark zero-cost windows
    for k = 1:numel(window_times_hr)
        xline(ax1, window_times_hr(k), '--', 'Color', greenC, 'LineWidth', 1.2, ...
              'Alpha', 0.8);
    end

    ylabel(ax1, '\DeltaV dog-leg (m/s)', 'Color', txtCol);
    title(ax1, title_str, 'Color', txtCol, 'FontSize', 10);
    xlim(ax1, [0 t_hr(end)]);
    ylim(ax1, [0 max_dv_m_s * 1.15 + 1]);

    if ~isempty(window_times_hr)
        text(ax1, window_times_hr(1), max_dv_m_s * 0.95, 'Natural windows', ...
             'Color', greenC, 'FontSize', 8);
    end

    %% Bottom panel: RAAN achieved vs time
    ax2 = subplot(2, 1, 2);
    set(ax2, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
             'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
    hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

    plot(ax2, t_hr, RAAN_asc,  '-',  'Color', cyanC, 'LineWidth', 1.5, 'DisplayName', 'Ascending');
    plot(ax2, t_hr, RAAN_desc, '-',  'Color', orngC, 'LineWidth', 1.5, 'DisplayName', 'Descending');
    yline(ax2, target_RAAN, '--', 'Color', [0.90 0.90 0.30], 'LineWidth', 1.5, ...
          'Label', sprintf('Target RAAN=%.1f°', target_RAAN), ...
          'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', ...
          'FontSize', 8);

    ylabel(ax2, 'Achieved RAAN (deg)', 'Color', txtCol);
    xlabel(ax2, 'Time from start (hours)', 'Color', txtCol);
    xlim(ax2, [0 t_hr(end)]);
    ylim(ax2, [0 360]);
    yticks(ax2, 0:60:360);

    leg = legend(ax2, 'Ascending', 'Descending', 'Location', 'northeast');
    set(leg, 'TextColor', txtCol, 'Color', bgCol + 0.05, 'EdgeColor', [0.4 0.4 0.4]);

    linkaxes([ax1, ax2], 'x');
end
end
