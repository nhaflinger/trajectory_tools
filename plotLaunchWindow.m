function fig = plotLaunchWindow(wins, launch_lat, launch_lon, target_inc, target_RAAN, start_jd, n_days)
%PLOTLAUNCHWINDOW  Visualize launch window analysis.
%
%   fig = plotLaunchWindow(wins, launch_lat, launch_lon, target_inc, target_RAAN, start_jd, n_days)
%
%   wins         - struct array returned by launchWindow()
%   launch_lat   - launch site latitude (deg)
%   launch_lon   - launch site longitude (deg)
%   target_inc   - target inclination (deg)
%   target_RAAN  - target RAAN at start_jd (deg)
%   start_jd     - start Julian Date (used for time axis)
%   n_days       - number of days in search window
%
%   Two-panel figure (dark theme):
%     Top:    LMST vs time — actual site LMST curve, required LMST lines, window markers
%     Bottom: Launch azimuth for each window — ascending (cyan) vs descending (orange)

%% ── Colours / style ──────────────────────────────────────────────────────────
bgCol  = [0.10 0.12 0.16];
axCol  = [0.14 0.18 0.26];
txtCol = [0.88 0.88 0.88];
gridC  = [0.25 0.28 0.35];
colAsc  = [0.20 0.85 0.95];   % cyan   — ascending windows
colDesc = [1.00 0.60 0.20];   % orange — descending windows
colReq  = [0.70 0.70 0.70];   % grey   — required LMST lines

%% ── Time vector ──────────────────────────────────────────────────────────────
t_hr = linspace(0, n_days*24, 1440*n_days + 1);    % 1-minute resolution, in hours
jd_vec = start_jd + t_hr / 24;

GMST_vec = mod(280.46061837 + 360.98564736629*(jd_vec - 2451545.0), 360);
LMST_vec = mod(GMST_vec + launch_lon, 360);

% Required LMST for ascending and descending windows
% Use the same geometry as launchWindow.m
phi = launch_lat;
cos_az = max(-1, min(1, cosd(target_inc) / cosd(phi)));
az_asc  = asind(cos_az);
az_desc = 180 - az_asc;

sin_uL = max(-1, min(1, sind(phi) / sind(target_inc)));
u_L = asind(sin_uL);
u_L_desc = 180 - u_L;

delta_lon_asc  = atan2d(tand(u_L)      * cosd(target_inc), 1);
delta_lon_desc = atan2d(tand(u_L_desc) * cosd(target_inc), 1);

LMST_req_asc  = mod(target_RAAN + delta_lon_asc,  360);
LMST_req_desc = mod(target_RAAN + delta_lon_desc, 360);

%% ── Figure ───────────────────────────────────────────────────────────────────
fig = figure('Color', bgCol, 'Position', [80 80 1200 700]);

titleStr = sprintf('Launch Windows  |  Site: lat=%.3f°, lon=%.3f°  |  Target: i=%.2f°, RAAN=%.2f°', ...
                   launch_lat, launch_lon, target_inc, target_RAAN);

%% ── Top panel: LMST vs Time ──────────────────────────────────────────────────
ax1 = subplot(2, 1, 1, 'Parent', fig);
set(ax1, 'Color', axCol, 'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', gridC, 'GridAlpha', 0.5, 'FontSize', 9.5);
hold(ax1, 'on');  grid(ax1, 'on');

% LMST advances as a sawtooth — plot in segments between wraps
% Break at 360 -> 0 transitions
lmst_plot = LMST_vec;
brk = find(diff(lmst_plot) < -180);
seg_starts = [1, brk+1];
seg_ends   = [brk, numel(lmst_plot)];

for s = 1:numel(seg_starts)
    idx = seg_starts(s):seg_ends(s);
    plot(ax1, t_hr(idx), lmst_plot(idx), '-', ...
        'Color', [0.30 0.65 0.98], 'LineWidth', 1.2, ...
        'HandleVisibility', iif(s==1,'on','off'), 'DisplayName', 'LMST at site');
end

% Required LMST lines (horizontal for no RAAN drift; drawn across full time axis)
yline(ax1, LMST_req_asc,  '--', 'Color', colAsc*0.85,  'LineWidth', 1.2, ...
    'DisplayName', sprintf('Req. LMST (ascending, %.1f deg)', LMST_req_asc));
yline(ax1, LMST_req_desc, '--', 'Color', colDesc*0.85, 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Req. LMST (descending, %.1f deg)', LMST_req_desc));

% Mark launch windows
for k = 1:numel(wins)
    w   = wins(k);
    tW  = (w.jd - start_jd) * 24;
    col = colAsc;
    lbl = 'off';
    if strcmpi(w.type, 'descending'), col = colDesc; end
    if k == 1, lbl = 'on'; end
    scatter(ax1, tW, w.LMST_deg, 70, 'o', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 0.7, ...
        'DisplayName', sprintf('Window %d (%s)', k, w.type), ...
        'HandleVisibility', lbl);
    % Label each window
    text(ax1, tW + 0.1, w.LMST_deg + 8, sprintf('%d', k), ...
        'Color', col, 'FontSize', 8, 'HorizontalAlignment', 'left');
end

xlabel(ax1, 'Hours from Search Start', 'Color', txtCol);
ylabel(ax1, 'LMST at Launch Site (deg)', 'Color', txtCol);
title(ax1, titleStr, 'Color', txtCol, 'FontSize', 10);
xlim(ax1, [0, n_days*24]);
ylim(ax1, [0, 360]);
yticks(ax1, 0:60:360);
legend(ax1, 'Location', 'northeastoutside', 'TextColor', txtCol, ...
    'Color', bgCol, 'EdgeColor', gridC, 'FontSize', 8);

%% ── Bottom panel: Azimuth bar chart ─────────────────────────────────────────
ax2 = subplot(2, 1, 2, 'Parent', fig);
set(ax2, 'Color', axCol, 'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', gridC, 'GridAlpha', 0.5, 'FontSize', 9.5);
hold(ax2, 'on');  grid(ax2, 'on');

n_wins = numel(wins);
if n_wins > 0
    win_times = ([wins.jd] - start_jd) * 24;
    win_az    = [wins.azimuth_deg];
    win_types = {wins.type};

    for k = 1:n_wins
        if strcmpi(win_types{k}, 'ascending')
            col = colAsc;
        else
            col = colDesc;
        end
        bar(ax2, win_times(k), win_az(k), 0.4, 'FaceColor', col, 'EdgeColor', 'none');
        text(ax2, win_times(k), win_az(k) + 1.5, sprintf('%d', k), ...
            'Color', txtCol, 'FontSize', 8, 'HorizontalAlignment', 'center');
    end

    % Legend proxies
    bar(ax2, NaN, NaN, 'FaceColor', colAsc,  'EdgeColor', 'none', 'DisplayName', 'Ascending');
    bar(ax2, NaN, NaN, 'FaceColor', colDesc, 'EdgeColor', 'none', 'DisplayName', 'Descending');

    ylim(ax2, [0, 200]);
    xlim(ax2, [0, n_days*24]);
else
    text(0.5, 0.5, 'No windows found', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', txtCol, 'FontSize', 12);
    xlim(ax2, [0, n_days*24]);  ylim(ax2, [0, 200]);
end

% Reference lines for north (0), east (90), south (180)
yline(ax2, 90,  ':', 'Color', txtCol*0.5, 'LineWidth', 0.8, 'HandleVisibility', 'off');
yline(ax2, 180, ':', 'Color', txtCol*0.5, 'LineWidth', 0.8, 'HandleVisibility', 'off');
text(ax2, n_days*24*0.01,  92, 'East (90°)',  'Color', txtCol*0.7, 'FontSize', 7.5);
text(ax2, n_days*24*0.01, 182, 'South (180°)', 'Color', txtCol*0.7, 'FontSize', 7.5);

xlabel(ax2, 'Hours from Search Start', 'Color', txtCol);
ylabel(ax2, 'Launch Azimuth (deg from North)', 'Color', txtCol);
title(ax2, 'Launch Azimuth per Window', 'Color', txtCol, 'FontSize', 10);
legend(ax2, 'Location', 'northeastoutside', 'TextColor', txtCol, ...
    'Color', bgCol, 'EdgeColor', gridC, 'FontSize', 8);

end

%% ── Local helpers ─────────────────────────────────────────────────────────────
function out = iif(cond, a, b)
    if cond, out = a; else, out = b; end
end
