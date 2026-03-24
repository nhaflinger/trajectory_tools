function fig = plotBetaAngle(result)
%PLOTBETAANGLE  Plot beta angle history and eclipse duration.
%
%   fig = plotBetaAngle(result)
%
%   Input:
%     result - struct returned by betaAngle()
%
%   Output:
%     fig    - figure handle
%
%   Two-panel dark-theme figure:
%     Top   : beta angle (deg) vs time (days), with eclipse boundary shaded
%     Bottom: eclipse duration per orbit (min) vs time (days)

bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];
accent = [0.30 0.75 0.93];   % cyan
warn   = [0.95 0.60 0.20];   % orange
shade  = [0.90 0.30 0.30];   % red for eclipse zone

orb = result.orb;
rho = result.rho_deg;

%% ── Title string ─────────────────────────────────────────────────────────────
alt_km = orb.a - 6378.1363;
title_str = sprintf('%s orbit | alt=%.0f km | i=%.1f deg | %.0f days', ...
    upper(orb.type), alt_km, orb.i, result.t_days(end));

%% ── Figure setup ─────────────────────────────────────────────────────────────
fig = figure('Color', bgCol, 'Position', [100 100 900 600]);

%% ── Top panel: beta angle ────────────────────────────────────────────────────
ax1 = subplot(2, 1, 1);
set(ax1, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
         'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

% Shade eclipse zone: |beta| < rho
t_max = result.t_days(end);
patch(ax1, [0 t_max t_max 0], [rho rho 90 90], shade, ...
      'FaceAlpha', 0.12, 'EdgeColor', 'none');
patch(ax1, [0 t_max t_max 0], [-rho -rho -90 -90], shade, ...
      'FaceAlpha', 0.12, 'EdgeColor', 'none');

% Eclipse boundary lines
plot(ax1, [0 t_max], [ rho  rho], '--', 'Color', shade, 'LineWidth', 1.2);
plot(ax1, [0 t_max], [-rho -rho], '--', 'Color', shade, 'LineWidth', 1.2);

% Beta angle
plot(ax1, result.t_days, result.beta_deg, '-', 'Color', accent, 'LineWidth', 1.5);

% Annotations
text(ax1, t_max * 0.02, rho + 2, sprintf('Eclipse zone (\\rho = %.1f°)', rho), ...
     'Color', shade, 'FontSize', 8);

ylabel(ax1, 'Beta angle (deg)', 'Color', txtCol);
title(ax1, title_str, 'Color', txtCol, 'FontSize', 11);
ylim(ax1, [-90 90]);
xlim(ax1, [0 t_max]);

% Stats annotation
stats_str = sprintf('max: %+.1f°  min: %+.1f°', ...
    result.beta_max_deg, result.beta_min_deg);
text(ax1, t_max * 0.02, 80, stats_str, 'Color', txtCol, 'FontSize', 8);

if result.eclipse_free
    text(ax1, t_max * 0.98, 80, 'ECLIPSE-FREE', ...
         'Color', [0.30 0.90 0.40], 'FontSize', 9, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right');
end

%% ── Bottom panel: eclipse duration ──────────────────────────────────────────
ax2 = subplot(2, 1, 2);
set(ax2, 'Color', bgCol, 'XColor', txtCol, 'YColor', txtCol, ...
         'GridColor', [0.35 0.35 0.40], 'GridAlpha', 0.4);
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

% Fill area under eclipse curve
patch(ax2, [result.t_days; flipud(result.t_days)], ...
           [result.t_eclipse_min; zeros(size(result.t_eclipse_min))], ...
      warn, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
plot(ax2, result.t_days, result.t_eclipse_min, '-', 'Color', warn, 'LineWidth', 1.5);

ylabel(ax2, 'Eclipse duration (min/orbit)', 'Color', txtCol);
xlabel(ax2, 'Time (days from epoch)', 'Color', txtCol);
xlim(ax2, [0 t_max]);
ylim(ax2, [0 max(result.t_eclipse_min) * 1.15 + 0.1]);

% Max eclipse annotation
[max_ecl, idx_max] = max(result.t_eclipse_min);
if max_ecl > 0
    text(ax2, result.t_days(idx_max), max_ecl * 1.05, ...
         sprintf('max: %.1f min', max_ecl), ...
         'Color', txtCol, 'FontSize', 8, 'HorizontalAlignment', 'center');
end

%% ── Link x-axes ──────────────────────────────────────────────────────────────
linkaxes([ax1, ax2], 'x');
end
