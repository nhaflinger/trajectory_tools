function fig = plotCoverage(cov, varargin)
%PLOTCOVERAGE  Plot coverage analysis results on a world map.
%
%   fig = plotCoverage(cov)
%   fig = plotCoverage(cov, Name, Value, ...)
%
%   Input:
%     cov - struct returned by coverageAnalysis()
%
%   Options (Name-Value pairs):
%     'Quantity'         - 'coverage' (default) or 'revisit'
%     'GroundStations'   - struct array with .lat .lon .name fields

p = inputParser;
addParameter(p, 'Quantity',       'coverage', @ischar);
addParameter(p, 'GroundStations', [],         @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts = p.Results;

%% ── Style constants (dark theme) ───────────────────────────────────────────
bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];
axCol  = [0.14 0.18 0.26];

%% ── Select quantity to plot ────────────────────────────────────────────────
quantity = lower(strtrim(opts.Quantity));
switch quantity
    case 'coverage'
        map_data    = cov.coverage_frac;
        cmap        = parula;
        clim_vals   = [0 1];
        cbar_label  = 'Coverage Fraction';
        hist_xlabel = 'Coverage Fraction';
        hist_data   = map_data(:);
        hist_xlim   = [0 1];
    case 'revisit'
        map_data    = cov.revisit_mean_hr;
        % Replace NaN with 0 for histogram but keep NaN on map
        hist_data   = map_data(~isnan(map_data));
        cmap        = flipud(hot);
        clim_vals   = [0, max(1, max(map_data(~isnan(map_data(:)))))];
        cbar_label  = 'Mean Revisit Time (hr)';
        hist_xlabel = 'Mean Revisit Time (hr)';
        hist_xlim   = [0, clim_vals(2)];
    otherwise
        error('plotCoverage: unknown Quantity ''%s''. Valid: coverage, revisit', quantity);
end

%% ── Coastlines ─────────────────────────────────────────────────────────────
[coastlon, coastlat] = loadCoastlines();

%% ── Figure layout ──────────────────────────────────────────────────────────
fig = figure('Color', bgCol, 'Position', [60 60 1300 780]);

%% ── Top panel: map ─────────────────────────────────────────────────────────
ax_map = subplot(2, 1, 1, 'Parent', fig);
set(ax_map, 'Color', axCol, ...
    'XColor', txtCol * 0.7, 'YColor', txtCol * 0.7, ...
    'GridColor', [0.25 0.28 0.35], 'GridAlpha', 0.4);
hold(ax_map, 'on');
grid(ax_map, 'on');

% pcolor map (lat x lon grid)
[LON_g, LAT_g] = meshgrid(cov.lon_vec, cov.lat_vec);
h_pc = pcolor(ax_map, LON_g, LAT_g, map_data);
set(h_pc, 'EdgeColor', 'none');
shading(ax_map, 'flat');

colormap(ax_map, cmap);
clim(ax_map, clim_vals);

cb = colorbar(ax_map);
cb.Color = txtCol;
cb.Label.String = cbar_label;
cb.Label.Color  = txtCol;
cb.Label.FontSize = 10;

% Coastlines
if ~isempty(coastlon)
    plot(ax_map, coastlon, coastlat, '-', ...
        'Color', [0.65 0.72 0.80], 'LineWidth', 0.7, 'HandleVisibility', 'off');
end

% Ground stations
if ~isempty(opts.GroundStations)
    for gs = opts.GroundStations(:)'
        scatter(ax_map, gs.lon, gs.lat, 80, 's', ...
            'MarkerFaceColor', [1.0 0.50 0.10], ...
            'MarkerEdgeColor', 'none', ...
            'DisplayName', gs.name);
        text(ax_map, gs.lon + 2, gs.lat + 3, gs.name, ...
            'Color', [1.0 0.72 0.35], 'FontSize', 8);
    end
end

xlim(ax_map, [-180 180]);
ylim(ax_map, [-90 90]);
xticks(ax_map, -180:30:180);
yticks(ax_map, -90:30:90);
xlabel(ax_map, 'Longitude (deg)', 'Color', txtCol, 'FontSize', 10);
ylabel(ax_map, 'Latitude (deg)',  'Color', txtCol, 'FontSize', 10);

orb = cov.orb;
title_str = sprintf('%s  |  alt: %.0f–%.0f km  |  i = %.2f°  |  e = %.4f  |  Duration: %.1f hr', ...
    upper(orb.type), orb.alt_peri, orb.alt_apo, orb.i, orb.e, cov.duration_hr);
title(ax_map, title_str, 'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');

%% ── Bottom panel: histogram ─────────────────────────────────────────────────
ax_hist = subplot(2, 1, 2, 'Parent', fig);
set(ax_hist, 'Color', axCol, ...
    'XColor', txtCol * 0.7, 'YColor', txtCol * 0.7, ...
    'GridColor', [0.25 0.28 0.35], 'GridAlpha', 0.4);
hold(ax_hist, 'on');
grid(ax_hist, 'on');

if ~isempty(hist_data)
    n_bins  = min(50, max(10, round(sqrt(numel(hist_data)))));
    h_hist  = histogram(ax_hist, hist_data, n_bins, ...
        'FaceColor', [0.25 0.72 0.98], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.85);

    if strcmp(quantity, 'coverage')
        xline(ax_hist, mean(hist_data, 'omitnan'), '--', ...
            'Color', [1.0 0.85 0.25], 'LineWidth', 1.5, ...
            'Label', sprintf('Mean = %.3f', mean(hist_data, 'omitnan')), ...
            'LabelColor', [1.0 0.85 0.25], 'FontSize', 9);
    else
        if ~isempty(hist_data)
            xline(ax_hist, mean(hist_data, 'omitnan'), '--', ...
                'Color', [1.0 0.85 0.25], 'LineWidth', 1.5, ...
                'Label', sprintf('Mean = %.2f hr', mean(hist_data, 'omitnan')), ...
                'LabelColor', [1.0 0.85 0.25], 'FontSize', 9);
        end
    end
end

xlim(ax_hist, hist_xlim);
xlabel(ax_hist, hist_xlabel, 'Color', txtCol, 'FontSize', 10);
ylabel(ax_hist, 'Number of Grid Points', 'Color', txtCol, 'FontSize', 10);
title(ax_hist, ['Distribution of ' cbar_label], 'Color', txtCol, 'FontSize', 10);
end

%% ── Local: coastline loader (copy from plotGroundTrack) ─────────────────────
function [lon, lat] = loadCoastlines()
%LOADCOASTLINES  Load coastline data without requiring the Mapping Toolbox.
    lon = [];  lat = [];
    cachePath = fullfile(fileparts(mfilename('fullpath')), 'coastlines_cache.mat');
    if isfile(cachePath)
        s = load(cachePath, 'coastlon', 'coastlat');
        lon = s.coastlon;  lat = s.coastlat;
        return;
    end
    try
        s = load('coastlines');   % Mapping Toolbox
        lon = s.coastlon;  lat = s.coastlat;
        return;
    catch
    end
    fprintf('plotCoverage: coastlines_cache.mat not found. Downloading now...\n');
    try
        downloadCoastlines();
        s = load(cachePath, 'coastlon', 'coastlat');
        lon = s.coastlon;  lat = s.coastlat;
    catch ME
        warning('plotCoverage:noCoastlines', ...
            'Could not load or download coastlines: %s\nPlotting without coastlines.', ME.message);
    end
end
