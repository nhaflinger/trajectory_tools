function [fig, summary] = plotAccessWindowGantt(orb, ground_stations, duration_hr, varargin)
%PLOTACCESSWINDOWGANTT  Gantt-style access window timeline for multiple ground stations.
%
%   fig = plotAccessWindowGantt(orb, ground_stations, duration_hr)
%   fig = plotAccessWindowGantt(orb, ground_stations, duration_hr, Name, Value, ...)
%   [fig, summary] = plotAccessWindowGantt(...)
%
%   Inputs:
%     orb             - earthOrbit struct
%     ground_stations - struct array, each with fields:
%                         .lat   geodetic latitude (deg)
%                         .lon   longitude (deg)
%                         .alt   altitude (km), optional, default 0
%                         .name  label string
%     duration_hr     - analysis duration (hours)
%
%   Options (Name-Value pairs):
%     'MinElevation'     - minimum elevation for contact (deg), default 10
%     'StepSize'         - propagation step size (s), default 30
%     'ColorByElevation' - color bars by max elevation (logical), default true
%     'ShowMaxEl'        - annotate bars with max elevation text, default true
%     'TimeFormat'       - 'hours' or 'minutes' for x-axis, default 'hours'
%
%   Output:
%     fig     - figure handle
%     summary - struct array (one per station) with fields:
%                 name, n_passes, total_contact_min, mean_duration_min,
%                 max_elevation_deg, windows

%% ── Parse inputs ────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'MinElevation',     10,     @isnumeric);
addParameter(p, 'StepSize',         30,     @isnumeric);
addParameter(p, 'ColorByElevation', true,   @islogical);
addParameter(p, 'ShowMaxEl',        true,   @islogical);
addParameter(p, 'TimeFormat',       'hours',@ischar);
parse(p, varargin{:});
opts = p.Results;

minEl      = opts.MinElevation;
stepSize   = opts.StepSize;
colorByEl  = opts.ColorByElevation;
showMaxEl  = opts.ShowMaxEl;
timeFmt    = lower(strtrim(opts.TimeFormat));

duration_s = duration_hr * 3600;
n_stations = numel(ground_stations);

%% ── Propagate orbit once ────────────────────────────────────────────────────
res = propagateOrbit(orb, duration_s, 'Method', 'j2', 'StepSize', stepSize);

t_vec  = res.t;           % Nx1  (seconds from epoch)
r_eci  = res.r_eci;       % 3xN
jd_vec = orb.epoch_jd + t_vec / 86400;   % 1xN Julian dates

%% ── Convert time axis ───────────────────────────────────────────────────────
if strcmp(timeFmt, 'minutes')
    t_axis  = t_vec / 60;
    t_label = 'Time from Epoch (minutes)';
    t_scale = 1/60;
else
    t_axis  = t_vec / 3600;
    t_label = 'Time from Epoch (hours)';
    t_scale = 1/3600;
end

%% ── Compute passes for each ground station ──────────────────────────────────
summary(n_stations) = struct('name', '', 'n_passes', 0, ...
    'total_contact_min', 0, 'mean_duration_min', 0, ...
    'max_elevation_deg', 0, 'windows', []);

all_passes = cell(n_stations, 1);

for si = 1:n_stations
    gs = ground_stations(si);

    % Default altitude to 0 if missing
    if ~isfield(gs, 'alt') || isempty(gs.alt)
        obs_alt = 0;
    else
        obs_alt = gs.alt;
    end

    % Compute elevation time series using topocentricAzEl
    [~, el, ~] = topocentricAzEl(gs.lat, gs.lon, obs_alt, r_eci, jd_vec');

    % Edge-detection: find rises and falls
    in_view     = el > minEl;
    transitions = diff([0, in_view, 0]);
    rise_idx    = find(transitions ==  1);   % start of each pass
    fall_idx    = find(transitions == -1) - 1; % end of each pass (inclusive)

    n_passes = min(numel(rise_idx), numel(fall_idx));

    % Build windows struct array
    empty_win = struct('start_s', {}, 'stop_s', {}, 'duration_s', {}, ...
                       'max_elev_deg', {});
    wins = empty_win;

    for w = 1:n_passes
        i0 = rise_idx(w);
        i1 = fall_idx(w);
        if i0 > i1
            continue;
        end
        t0 = t_vec(i0);
        t1 = t_vec(i1);
        max_el_w = max(el(i0:i1));

        wins(end+1).start_s      = t0;   %#ok<AGROW>
        wins(end).stop_s         = t1;
        wins(end).duration_s     = t1 - t0;
        wins(end).max_elev_deg   = max_el_w;
    end

    all_passes{si} = wins;

    % Compute summary statistics
    n_w = numel(wins);
    if n_w > 0
        total_s  = sum([wins.duration_s]);
        mean_s   = total_s / n_w;
        max_el_s = max([wins.max_elev_deg]);
    else
        total_s  = 0;
        mean_s   = 0;
        max_el_s = 0;
    end

    summary(si).name               = gs.name;
    summary(si).n_passes           = n_w;
    summary(si).total_contact_min  = total_s / 60;
    summary(si).mean_duration_min  = mean_s  / 60;
    summary(si).max_elevation_deg  = max_el_s;
    summary(si).windows            = wins;
end

%% ── Print contact statistics table ─────────────────────────────────────────
fprintf('\n=== Access Window Summary ===\n');
hdr = sprintf('%-14s | %6s | %20s | %20s | %12s', ...
    'Station', 'Passes', 'Total Contact (min)', 'Mean Duration (min)', 'Max El (deg)');
fprintf('%s\n', hdr);
fprintf('%s\n', repmat('-', 1, numel(hdr)));
for si = 1:n_stations
    s = summary(si);
    fprintf('%-14s | %6d | %20.1f | %20.1f | %12.1f\n', ...
        s.name, s.n_passes, s.total_contact_min, s.mean_duration_min, s.max_elevation_deg);
end
fprintf('\n');

%% ── Build figure ────────────────────────────────────────────────────────────
bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];
axCol  = [0.14 0.18 0.26];

fig = figure('Color', bgCol, 'Position', [60 60 1200 max(380, 80 + n_stations*70)]);
ax  = axes('Parent', fig, 'Color', axCol, ...
    'XColor', txtCol * 0.85, 'YColor', txtCol * 0.85, ...
    'GridColor', [0.22 0.26 0.34], 'GridAlpha', 0.55, ...
    'TickDir', 'out', 'Box', 'on');
hold(ax, 'on');
grid(ax, 'on');

% ── Colormap for elevation ────────────────────────────────────────────────
if colorByEl
    cmap  = hot(256);
    cmap  = cmap(30:end, :);   % trim very dark end for readability
    n_cm  = size(cmap, 1);
    el_lo = minEl;
    el_hi = 90;
else
    station_cmap = lines(n_stations);
end

% Minimum bar width threshold for annotation
bar_thresh_s = duration_s / 50;

% ── Draw vertical dashed lines at each orbital period ────────────────────
T_orb = orb.period;   % seconds
period_times = T_orb : T_orb : duration_s;
for kp = 1:numel(period_times)
    xp = period_times(kp) * t_scale;
    line(ax, [xp xp], [0.4 n_stations + 0.6], ...
        'Color', [0.35 0.40 0.50], 'LineStyle', '--', 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
end

% ── Draw horizontal grid lines between stations ───────────────────────────
for si = 1:n_stations - 1
    y_line = si + 0.5;
    line(ax, [0, duration_s * t_scale], [y_line, y_line], ...
        'Color', [0.22 0.26 0.34], 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

% ── Draw pass bars ────────────────────────────────────────────────────────
bar_height = 0.62;

for si = 1:n_stations
    wins = all_passes{si};
    y_ctr = si;   % station row (y = 1 .. n_stations)

    for w = 1:numel(wins)
        x0 = wins(w).start_s * t_scale;
        x1 = wins(w).stop_s  * t_scale;
        dx = x1 - x0;
        el_w = wins(w).max_elev_deg;

        if colorByEl
            % Map max elevation to colormap index
            frac  = (el_w - el_lo) / (el_hi - el_lo);
            frac  = max(0, min(1, frac));
            ci    = max(1, round(1 + frac * (n_cm - 1)));
            bar_c = cmap(ci, :);
        else
            bar_c = station_cmap(si, :);
        end

        % Draw filled bar
        xv = [x0, x1, x1, x0];
        yv = [y_ctr - bar_height/2, y_ctr - bar_height/2, ...
              y_ctr + bar_height/2, y_ctr + bar_height/2];
        fill(ax, xv, yv, bar_c, 'EdgeColor', bar_c * 0.75, ...
            'LineWidth', 0.5, 'HandleVisibility', 'off');

        % Annotate with max elevation if bar is wide enough
        if showMaxEl && wins(w).duration_s > bar_thresh_s
            text(ax, x0 + dx/2, y_ctr, sprintf('%d\x00B0', round(el_w)), ...
                'Color', [0.08 0.06 0.04], 'FontSize', 7.5, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Clipping', 'on');
        end
    end
end

% ── Colorbar ─────────────────────────────────────────────────────────────
if colorByEl
    colormap(ax, cmap);
    cb = colorbar(ax, 'Location', 'eastoutside');
    caxis(ax, [el_lo, el_hi]);
    cb.Label.String    = 'Max Elevation (deg)';
    cb.Label.Color     = txtCol;
    cb.Color           = txtCol * 0.85;
    cb.FontSize        = 9;
end

% ── Summary text box ──────────────────────────────────────────────────────
% Build multi-line string of total contact per station
sumLines = cell(n_stations, 1);
for si = 1:n_stations
    s = summary(si);
    sumLines{si} = sprintf('%s: %d passes, %.0f min', ...
        s.name, s.n_passes, s.total_contact_min);
end
sumTxt = strjoin(sumLines, '   |   ');
text(ax, 0.01, 0.013, sumTxt, ...
    'Units', 'normalized', ...
    'Color', txtCol * 0.80, ...
    'FontSize', 7.5, ...
    'Interpreter', 'none', ...
    'BackgroundColor', [0.08 0.10 0.15], ...
    'EdgeColor', [0.30 0.35 0.45], ...
    'Margin', 3);

% ── Axis formatting ───────────────────────────────────────────────────────
xlim(ax, [0, duration_hr * (strcmp(timeFmt,'hours') + ~strcmp(timeFmt,'hours') * 60)]);
ylim(ax, [0.4, n_stations + 0.6]);

yticks(ax, 1:n_stations);
yticklabels(ax, {ground_stations.name});
set(ax, 'YDir', 'normal', 'FontSize', 9.5);

xlabel(ax, t_label, 'Color', txtCol, 'FontSize', 10);
ylabel(ax, 'Ground Station', 'Color', txtCol, 'FontSize', 10);

% Title
orb_type_str = '';
if isfield(orb, 'type')
    orb_type_str = [upper(orb.type(1)) lower(orb.type(2:end)) ' '];
end
alt_str = '';
if isfield(orb, 'a')
    R_E = 6378.1363;
    alt_km = orb.a - R_E;
    alt_str = sprintf('%.0f km', alt_km);
end
inc_str = '';
if isfield(orb, 'i')
    inc_str = sprintf('%.1f', orb.i);
end

title_str = sprintf('Access Windows: %s%s | i=%s° | %.0f hr | MinEl=%d°', ...
    orb_type_str, alt_str, inc_str, duration_hr, minEl);
title(ax, title_str, 'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');

% ── Period label on vertical lines ───────────────────────────────────────
if ~isempty(period_times)
    first_xp = period_times(1) * t_scale;
    text(ax, first_xp, n_stations + 0.55, 'T_{orb}', ...
        'Color', [0.45 0.52 0.62], 'FontSize', 7, ...
        'HorizontalAlignment', 'center', 'Clipping', 'on');
end

hold(ax, 'off');

end
