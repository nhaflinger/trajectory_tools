function [result, varargout] = linkBudget(varargin)
%LINKBUDGET  Two-way link budget analysis: static snapshot or pass time-series.
%
%   Static mode (snapshot at a given slant range):
%     result = linkBudget(Name, Value, ...)
%     [result, fig] = linkBudget(..., 'Plot', true)
%
%   Pass mode (full time series over 3 orbital periods):
%     result = linkBudget(orb, gs, Name, Value, ...)
%     [result, fig] = linkBudget(orb, gs, ..., 'Plot', true)
%
%   Pass mode is triggered when the first argument is an orbit struct (has
%   field 'a') and the second argument is a ground station struct (has field
%   'lat').
%
%   Static mode Name-Value parameters:
%     'Freq_GHz'      - Carrier frequency (GHz)            [default: 2.0]
%     'P_tx_dBW'      - Transmit power (dBW)               [default: 0]
%     'G_tx_dBi'      - Transmit antenna gain (dBi)        [default: 6]
%     'G_rx_dBi'      - Receive antenna gain (dBi)         [default: 45]
%     'T_sys_K'       - System noise temperature (K)       [default: 135]
%     'DataRate_bps'  - Data rate (bps)                    [default: 1e6]
%     'ReqEbN0_dB'    - Required Eb/N0 for closure (dB)   [default: 10]
%     'L_atm_dB'      - Atmospheric + rain loss (dB)       [default: 0.3]
%     'L_point_dB'    - Pointing loss (dB)                 [default: 0]
%     'Range_km'      - Slant range (km)                   [default: 1000]
%     'Plot'          - Generate figure                    [default: false]
%
%   Pass mode additional Name-Value parameters:
%     'MinElevation'  - Minimum elevation angle (deg)      [default: 5]
%
%   Ground station struct fields (gs):
%     gs.lat          - Latitude (deg)
%     gs.lon          - Longitude (deg)
%     gs.alt          - Altitude (km)        [default: 0]
%     gs.name         - Label string         [default: 'GS']
%     gs.Freq_GHz, gs.P_tx_dBW, etc. — optional overrides for link params
%
%   Static output struct fields:
%     freq_GHz, P_tx_dBW, G_tx_dBi, G_rx_dBi, T_sys_K, DataRate_bps,
%     ReqEbN0_dB, L_atm_dB, L_point_dB, EIRP_dBW, lambda_m, FSPL_dB,
%     P_rx_dBW, N0_dBW_Hz, C_N0_dBHz, Eb_N0_dB, LinkMargin_dB,
%     range_km, max_range_km
%
%   Pass output struct fields:
%     t_s, jd, az_deg, el_deg, range_km, range_rate_km_s, doppler_Hz,
%     FSPL_dB, P_rx_dBW, C_N0_dBHz, Eb_N0_dB, LinkMargin_dB,
%     passes (struct array), static_budget, gs, orb

%% ── Detect mode ─────────────────────────────────────────────────────────────
passMode = false;
if nargin >= 2 && isstruct(varargin{1}) && isfield(varargin{1}, 'a') && ...
        isstruct(varargin{2}) && isfield(varargin{2}, 'lat')
    passMode = true;
    orb    = varargin{1};
    gs     = varargin{2};
    nvArgs = varargin(3:end);
else
    nvArgs = varargin;
end

%% ── Parse Name-Value arguments ──────────────────────────────────────────────
p = inputParser;
addParameter(p, 'Freq_GHz',     2.0,    @isnumeric);
addParameter(p, 'P_tx_dBW',     0,      @isnumeric);
addParameter(p, 'G_tx_dBi',     6,      @isnumeric);
addParameter(p, 'G_rx_dBi',     45,     @isnumeric);
addParameter(p, 'T_sys_K',      135,    @isnumeric);
addParameter(p, 'DataRate_bps', 1e6,    @isnumeric);
addParameter(p, 'ReqEbN0_dB',   10,     @isnumeric);
addParameter(p, 'L_atm_dB',     0.3,    @isnumeric);
addParameter(p, 'L_point_dB',   0,      @isnumeric);
addParameter(p, 'Range_km',     1000,   @isnumeric);
addParameter(p, 'Plot',         false,  @(x) islogical(x) || isnumeric(x));
addParameter(p, 'MinElevation', 5,      @isnumeric);
addParameter(p, 'NumOrbits',   16,     @isnumeric);   % pass mode: orbits to propagate (~24 hr for LEO)
parse(p, nvArgs{:});
opts = p.Results;

% Override with ground station fields if present (pass mode)
if passMode
    fields = {'Freq_GHz','P_tx_dBW','G_tx_dBi','G_rx_dBi','T_sys_K', ...
              'DataRate_bps','ReqEbN0_dB','L_atm_dB','L_point_dB'};
    for fi = 1:numel(fields)
        fn = fields{fi};
        if isfield(gs, fn)
            opts.(fn) = gs.(fn);
        end
    end
    if ~isfield(gs, 'alt'),  gs.alt  = 0;    end
    if ~isfield(gs, 'name'), gs.name = 'GS'; end
end

doPlot = logical(opts.Plot);

if ~passMode
    %% ── STATIC MODE ─────────────────────────────────────────────────────────
    result = computeStatic(opts);
    printStaticTable(result);

    if doPlot
        fig = plotStaticBudget(result, opts);
        if nargout > 1
            varargout{1} = fig;
        end
    elseif nargout > 1
        varargout{1} = [];
    end
    return;
end

%% ── PASS MODE ───────────────────────────────────────────────────────────────
c_ms  = 2.998e8;   % speed of light (m/s)

% 1. Propagate for NumOrbits periods (default 16 ≈ 24 hr for LEO)
dur_s = opts.NumOrbits * orb.period;
traj  = propagateOrbit(orb, dur_s, 'Method', 'j2', 'StepSize', 30);

t_s   = traj.t;          % Nx1
r_eci = traj.r_eci;      % 3xN
N_t   = numel(t_s);

% 2. JD time vector
jd_vec = orb.epoch_jd + t_s / 86400;   % Nx1

% 3. Topocentric az/el/range at each step
[az_deg, el_deg, range_km] = topocentricAzEl( ...
    gs.lat, gs.lon, gs.alt, r_eci, jd_vec(:)');

az_deg   = az_deg(:);    % Nx1
el_deg   = el_deg(:);    % Nx1
range_km = range_km(:);  % Nx1

% 4. Range rate (km/s) — finite difference, forward, pad last element
dr              = diff(range_km);
dt_v            = diff(t_s);
rr              = dr ./ dt_v;
range_rate_km_s = [rr; rr(end)];   % Nx1

% 5. Doppler shift (Hz)
freq_Hz    = opts.Freq_GHz * 1e9;
doppler_Hz = -freq_Hz * range_rate_km_s * 1e3 / c_ms;   % Nx1

% 6. Per-timestep link budget quantities
lambda_m   = c_ms / freq_Hz;
EIRP_dBW   = opts.P_tx_dBW + opts.G_tx_dBi;
N0_dBW_Hz  = -228.6 + 10*log10(opts.T_sys_K);

range_m_all  = range_km * 1e3;
FSPL_dB_all  = 20*log10(4*pi .* range_m_all / lambda_m);
P_rx_dBW_all = EIRP_dBW - FSPL_dB_all - opts.L_atm_dB - opts.L_point_dB + opts.G_rx_dBi;
C_N0_all     = P_rx_dBW_all - N0_dBW_Hz;
Eb_N0_all    = C_N0_all - 10*log10(opts.DataRate_bps);
LM_all       = Eb_N0_all - opts.ReqEbN0_dB;

% 7. Mask below-horizon values to NaN so plots only show meaningful contact data
minEl = opts.MinElevation;
below_horizon = el_deg <= minEl;
FSPL_dB_all(below_horizon)  = NaN;
P_rx_dBW_all(below_horizon) = NaN;
C_N0_all(below_horizon)     = NaN;
Eb_N0_all(below_horizon)    = NaN;
LM_all(below_horizon)       = NaN;
doppler_Hz(below_horizon)   = NaN;

% 8. Detect passes (el > MinElevation)
above  = el_deg > minEl;

rising  = find(~above(1:end-1) &  above(2:end)) + 1;
falling = find( above(1:end-1) & ~above(2:end));

rising  = rising(:);
falling = falling(:);

if above(1),   rising  = [1;   rising];  end
if above(end), falling = [falling; N_t]; end

n_passes = min(numel(rising), numel(falling));

passes = struct('start_s', {}, 'stop_s', {}, 'duration_s', {}, ...
                'max_el_deg', {}, 'max_LinkMargin_dB', {}, ...
                'min_range_km', {}, 'contact_jd_start', {}, 'contact_jd_stop', {});

for w = 1:n_passes
    i0 = rising(w);
    i1 = falling(w);
    if i0 > i1, continue; end
    passes(end+1).start_s          = t_s(i0);       %#ok<AGROW>
    passes(end).stop_s             = t_s(i1);
    passes(end).duration_s         = t_s(i1) - t_s(i0);
    passes(end).max_el_deg         = max(el_deg(i0:i1));
    passes(end).max_LinkMargin_dB  = max(LM_all(i0:i1));
    passes(end).min_range_km       = min(range_km(i0:i1));
    passes(end).contact_jd_start   = jd_vec(i0);
    passes(end).contact_jd_stop    = jd_vec(i1);
end

% Static budget at minimum range of first pass (representative snapshot)
if ~isempty(passes)
    ref_range = passes(1).min_range_km;
else
    ref_range = opts.Range_km;
end
opts_static      = opts;
opts_static.Range_km = ref_range;
static_bud       = computeStatic(opts_static);

result = struct( ...
    't_s',             t_s,            ...
    'jd',              jd_vec,         ...
    'az_deg',          az_deg,         ...
    'el_deg',          el_deg,         ...
    'range_km',        range_km,       ...
    'range_rate_km_s', range_rate_km_s, ...
    'doppler_Hz',      doppler_Hz,     ...
    'FSPL_dB',         FSPL_dB_all,    ...
    'P_rx_dBW',        P_rx_dBW_all,   ...
    'C_N0_dBHz',       C_N0_all,       ...
    'Eb_N0_dB',        Eb_N0_all,      ...
    'LinkMargin_dB',   LM_all,         ...
    'passes',          passes,         ...
    'static_budget',   static_bud,     ...
    'gs',              gs,             ...
    'orb',             orb);

printPassTable(result, passes);

if doPlot
    fig = plotPassBudget(result, opts, orb, gs);
    if nargout > 1
        varargout{1} = fig;
    end
elseif nargout > 1
    varargout{1} = [];
end

end  % linkBudget


%% ════════════════════════════════════════════════════════════════════════════
%% ── Local: compute static link budget from opts struct ──────────────────────
function s = computeStatic(opts)
c_ms   = 2.998e8;
k_B_dB = -228.6;

freq_Hz  = opts.Freq_GHz * 1e9;
lambda_m = c_ms / freq_Hz;
range_m  = opts.Range_km * 1e3;

EIRP_dBW   = opts.P_tx_dBW + opts.G_tx_dBi;
FSPL_dB    = 20*log10(4*pi * range_m / lambda_m);
P_rx_dBW   = EIRP_dBW - FSPL_dB - opts.L_atm_dB - opts.L_point_dB + opts.G_rx_dBi;
N0_dBW_Hz  = k_B_dB + 10*log10(opts.T_sys_K);
C_N0_dBHz  = P_rx_dBW - N0_dBW_Hz;
Eb_N0_dB   = C_N0_dBHz - 10*log10(opts.DataRate_bps);
LM_dB      = Eb_N0_dB - opts.ReqEbN0_dB;

% Max closure range: solve LM = 0 analytically
% FSPL_budget = EIRP - L_atm - L_point + G_rx - N0 - 10log10(DR) - ReqEbN0
% Then range_max = (lambda/(4*pi)) * 10^(FSPL_budget/20)
FSPL_budget  = EIRP_dBW - opts.L_atm_dB - opts.L_point_dB + opts.G_rx_dBi ...
               - N0_dBW_Hz - 10*log10(opts.DataRate_bps) - opts.ReqEbN0_dB;
range_max_km = ((lambda_m / (4*pi)) * 10^(FSPL_budget / 20)) / 1e3;

s = struct( ...
    'freq_GHz',      opts.Freq_GHz,      ...
    'P_tx_dBW',      opts.P_tx_dBW,      ...
    'G_tx_dBi',      opts.G_tx_dBi,      ...
    'G_rx_dBi',      opts.G_rx_dBi,      ...
    'T_sys_K',       opts.T_sys_K,       ...
    'DataRate_bps',  opts.DataRate_bps,  ...
    'ReqEbN0_dB',    opts.ReqEbN0_dB,    ...
    'L_atm_dB',      opts.L_atm_dB,      ...
    'L_point_dB',    opts.L_point_dB,    ...
    'EIRP_dBW',      EIRP_dBW,           ...
    'lambda_m',      lambda_m,           ...
    'FSPL_dB',       FSPL_dB,            ...
    'P_rx_dBW',      P_rx_dBW,           ...
    'N0_dBW_Hz',     N0_dBW_Hz,          ...
    'C_N0_dBHz',     C_N0_dBHz,          ...
    'Eb_N0_dB',      Eb_N0_dB,           ...
    'LinkMargin_dB', LM_dB,              ...
    'range_km',      opts.Range_km,      ...
    'max_range_km',  range_max_km);
end

%% ── Local: print static table ───────────────────────────────────────────────
function printStaticTable(s)
fprintf('\n');
fprintf('+---------------------------------------------------------+\n');
fprintf('|            LINK BUDGET -- STATIC ANALYSIS               |\n');
fprintf('+---------------------------------------------------------+\n');
fprintf('| Frequency            %10.3f GHz                    |\n', s.freq_GHz);
fprintf('| Wavelength           %10.4f m                      |\n', s.lambda_m);
fprintf('| Slant Range          %10.1f km                     |\n', s.range_km);
fprintf('+---------------------------------------------------------+\n');
fprintf('| Transmit Power       %10.2f dBW                    |\n', s.P_tx_dBW);
fprintf('| Tx Antenna Gain      %10.2f dBi                    |\n', s.G_tx_dBi);
fprintf('| EIRP                 %10.2f dBW                    |\n', s.EIRP_dBW);
fprintf('+---------------------------------------------------------+\n');
fprintf('| Free-Space Path Loss %10.2f dB                     |\n', s.FSPL_dB);
fprintf('| Atmospheric Loss     %10.2f dB                     |\n', s.L_atm_dB);
fprintf('| Pointing Loss        %10.2f dB                     |\n', s.L_point_dB);
fprintf('+---------------------------------------------------------+\n');
fprintf('| Rx Antenna Gain      %10.2f dBi                    |\n', s.G_rx_dBi);
fprintf('| Received Power       %10.2f dBW                    |\n', s.P_rx_dBW);
fprintf('+---------------------------------------------------------+\n');
fprintf('| System Noise Temp    %10.1f K                      |\n', s.T_sys_K);
fprintf('| Noise PSD (N0)       %10.2f dBW/Hz                 |\n', s.N0_dBW_Hz);
fprintf('| C/N0                 %10.2f dBHz                   |\n', s.C_N0_dBHz);
fprintf('+---------------------------------------------------------+\n');
fprintf('| Data Rate            %10.0f bps                    |\n', s.DataRate_bps);
fprintf('| Eb/N0 (achieved)     %10.2f dB                     |\n', s.Eb_N0_dB);
fprintf('| Required Eb/N0       %10.2f dB                     |\n', s.ReqEbN0_dB);
fprintf('+---------------------------------------------------------+\n');
if s.LinkMargin_dB >= 0
    fprintf('| LINK MARGIN          %10.2f dB   LINK CLOSED        |\n', s.LinkMargin_dB);
else
    fprintf('| LINK MARGIN          %10.2f dB   LINK OPEN          |\n', s.LinkMargin_dB);
end
fprintf('| Max Closure Range    %10.1f km                     |\n', s.max_range_km);
fprintf('+---------------------------------------------------------+\n');
fprintf('\n');
end

%% ── Local: print pass summary table ────────────────────────────────────────
function printPassTable(result, passes)
fprintf('\n');
fprintf('PASS ANALYSIS SUMMARY -- %s -> %s\n', upper(result.orb.type), result.gs.name);
fprintf('%-6s  %-12s  %-10s  %-12s  %-14s\n', ...
    'Pass #', 'Duration (s)', 'MaxEl (deg)', 'MinRange (km)', 'LinkMargin (dB)');
fprintf('%s\n', repmat('-', 1, 62));
for w = 1:numel(passes)
    fprintf('  %-4d  %12.1f  %10.2f  %12.1f  %14.2f\n', w, ...
        passes(w).duration_s, passes(w).max_el_deg, ...
        passes(w).min_range_km, passes(w).max_LinkMargin_dB);
end
fprintf('\n');
end

%% ── Local: static mode plot ─────────────────────────────────────────────────
function fig = plotStaticBudget(s, opts)
bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];
axCol  = [0.14 0.18 0.26];
green  = [0.20 0.90 0.50];
red    = [0.95 0.35 0.35];

fig = figure('Color', bgCol, 'Position', [100 100 900 480]);
ax  = axes('Parent', fig, 'Color', axCol, ...
    'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', [0.25 0.28 0.35], 'GridAlpha', 0.5);
hold(ax, 'on');  grid(ax, 'on');

% Range sweep
r_vec   = logspace(log10(500), log10(50000), 500);
c_ms    = 2.998e8;
freq_Hz = opts.Freq_GHz * 1e9;
lambda  = c_ms / freq_Hz;
EIRP    = opts.P_tx_dBW + opts.G_tx_dBi;
N0      = -228.6 + 10*log10(opts.T_sys_K);
FSPL_v  = 20*log10(4*pi .* (r_vec * 1e3) / lambda);
Prx_v   = EIRP - FSPL_v - opts.L_atm_dB - opts.L_point_dB + opts.G_rx_dBi;
CN0_v   = Prx_v - N0;
EbN0_v  = CN0_v - 10*log10(opts.DataRate_bps);
LM_v    = EbN0_v - opts.ReqEbN0_dB;

% Shade positive / negative margin regions
pos_mask = LM_v >= 0;
neg_mask = ~pos_mask;

if any(pos_mask)
    fill(ax, [r_vec(pos_mask), fliplr(r_vec(pos_mask))], ...
         [LM_v(pos_mask), zeros(1, sum(pos_mask))], ...
         green, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
if any(neg_mask)
    fill(ax, [r_vec(neg_mask), fliplr(r_vec(neg_mask))], ...
         [LM_v(neg_mask), zeros(1, sum(neg_mask))], ...
         red, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Zero margin line
yline(ax, 0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0, ...
    'DisplayName', 'Zero margin');

% Link margin curve
plot(ax, r_vec, LM_v, '-', 'Color', [0.35 0.70 1.00], 'LineWidth', 2.0, ...
    'DisplayName', 'Link margin');

% Mark operating range
xline(ax, s.range_km, '--', 'Color', [0.95 0.85 0.25], 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
scatter(ax, s.range_km, s.LinkMargin_dB, 80, 'o', ...
    'MarkerFaceColor', [0.95 0.85 0.25], 'MarkerEdgeColor', 'none', ...
    'DisplayName', sprintf('Operating range %.0f km (LM = %.1f dB)', ...
    s.range_km, s.LinkMargin_dB));

% Mark zero-margin range
if s.max_range_km < 5e4 && s.max_range_km > 0
    xline(ax, s.max_range_km, ':', 'Color', red, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
    scatter(ax, s.max_range_km, 0, 80, '^', ...
        'MarkerFaceColor', red, 'MarkerEdgeColor', 'none', ...
        'DisplayName', sprintf('Max closure range %.0f km', s.max_range_km));
end

set(ax, 'XScale', 'log');
xlabel(ax, 'Slant Range (km)', 'Color', txtCol, 'FontSize', 11);
ylabel(ax, 'Link Margin (dB)', 'Color', txtCol, 'FontSize', 11);

GT_dBK = opts.G_rx_dBi - 10*log10(opts.T_sys_K);
title_str = sprintf('Link Budget  |  %.3f GHz  |  EIRP = %.1f dBW  |  G/T = %.1f dB/K', ...
    opts.Freq_GHz, EIRP, GT_dBK);
title(ax, title_str, 'Color', txtCol, 'FontSize', 11, 'FontWeight', 'bold');

legend(ax, 'Location', 'northeast', 'TextColor', txtCol, ...
    'Color', [0.08 0.10 0.14], 'EdgeColor', [0.30 0.32 0.38]);
end

%% ── Local: pass mode plot ───────────────────────────────────────────────────
function fig = plotPassBudget(result, opts, orb, gs)
bgCol  = [0.10 0.12 0.16];
txtCol = [0.88 0.88 0.88];
axCol  = [0.14 0.18 0.26];
gridC  = [0.25 0.28 0.35];

t_min  = result.t_s / 60;   % minutes from epoch
passes = result.passes;
n_pass = numel(passes);

if n_pass > 1
    pass_cmap = lines(n_pass);
else
    pass_cmap = [0.25 0.72 0.98];
end

dr_kbps    = opts.DataRate_bps / 1e3;
main_title = sprintf('Link Budget: %s  ->  %s  |  %.3f GHz  |  %.1f kbps', ...
    upper(orb.type), gs.name, opts.Freq_GHz, dr_kbps);

fig = figure('Color', bgCol, 'Position', [80 60 1300 820]);
annotation(fig, 'textbox', [0 0.96 1 0.04], 'String', main_title, ...
    'Color', txtCol, 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
    'BackgroundColor', bgCol);

%% Panel 1: Elevation vs time
ax1 = subplot(2, 2, 1, 'Parent', fig);
set(ax1, 'Color', axCol, 'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', gridC, 'GridAlpha', 0.5);
hold(ax1, 'on');  grid(ax1, 'on');

plot(ax1, t_min, result.el_deg, '-', 'Color', [0.35 0.70 1.00], 'LineWidth', 1.3, ...
    'HandleVisibility', 'off');
yline(ax1, opts.MinElevation, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.9, ...
    'HandleVisibility', 'off');

for w = 1:n_pass
    t0m = passes(w).start_s / 60;
    t1m = passes(w).stop_s  / 60;
    col = pass_cmap(min(w, size(pass_cmap,1)), :);
    patch(ax1, [t0m t1m t1m t0m], ...
          [opts.MinElevation opts.MinElevation 90 90], ...
          col, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    xline(ax1, t0m, '--', 'Color', col*0.8+0.05, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(ax1, t1m, '--', 'Color', col*0.8+0.05, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    % Pass label at peak elevation
    mask_w = result.t_s >= passes(w).start_s & result.t_s <= passes(w).stop_s;
    if any(mask_w)
        [max_el_w, idx_max] = max(result.el_deg(mask_w));
        t_max_w = t_min(mask_w);
        text(ax1, t_max_w(idx_max), max_el_w * 0.82, sprintf('P%d', w), ...
            'Color', col, 'FontSize', 8, 'HorizontalAlignment', 'center');
    end
end

xlabel(ax1, 'Time from Epoch (min)', 'Color', txtCol, 'FontSize', 10);
ylabel(ax1, 'Elevation (deg)',        'Color', txtCol, 'FontSize', 10);
title(ax1,  'Elevation Angle',        'Color', txtCol, 'FontSize', 10, 'FontWeight', 'bold');
ylim(ax1, [0 90]);

%% Panel 2: Link margin vs time
ax2 = subplot(2, 2, 2, 'Parent', fig);
set(ax2, 'Color', axCol, 'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', gridC, 'GridAlpha', 0.5);
hold(ax2, 'on');  grid(ax2, 'on');

plot(ax2, t_min, result.LinkMargin_dB, '-', 'Color', [0.35 0.70 1.00], 'LineWidth', 1.3, ...
    'HandleVisibility', 'off');
yline(ax2, 0, '--', 'Color', [0.95 0.35 0.35], 'LineWidth', 1.0, 'HandleVisibility', 'off');

for w = 1:n_pass
    t0m = passes(w).start_s / 60;
    t1m = passes(w).stop_s  / 60;
    col = pass_cmap(min(w, size(pass_cmap,1)), :);
    xline(ax2, t0m, '--', 'Color', col*0.8+0.05, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(ax2, t1m, '--', 'Color', col*0.8+0.05, 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

xlabel(ax2, 'Time from Epoch (min)', 'Color', txtCol, 'FontSize', 10);
ylabel(ax2, 'Link Margin (dB)',      'Color', txtCol, 'FontSize', 10);
title(ax2,  'Link Margin vs Time',   'Color', txtCol, 'FontSize', 10, 'FontWeight', 'bold');

%% Panel 3: Doppler shift during passes (kHz)
ax3 = subplot(2, 2, 3, 'Parent', fig);
set(ax3, 'Color', axCol, 'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', gridC, 'GridAlpha', 0.5);
hold(ax3, 'on');  grid(ax3, 'on');

for w = 1:n_pass
    col  = pass_cmap(min(w, size(pass_cmap,1)), :);
    mask = result.t_s >= passes(w).start_s & result.t_s <= passes(w).stop_s;
    if any(mask)
        plot(ax3, t_min(mask), result.doppler_Hz(mask) / 1e3, '-', ...
            'Color', col, 'LineWidth', 1.5, 'DisplayName', sprintf('Pass %d', w));
    end
end

if n_pass == 0
    text(ax3, 0.5, 0.5, 'No passes above minimum elevation', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', txtCol * 0.6, 'FontSize', 10);
end

yline(ax3, 0, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.9, 'HandleVisibility', 'off');
xlabel(ax3, 'Time from Epoch (min)', 'Color', txtCol, 'FontSize', 10);
ylabel(ax3, 'Doppler Shift (kHz)',   'Color', txtCol, 'FontSize', 10);
title(ax3,  'Doppler Shift During Passes', 'Color', txtCol, 'FontSize', 10, 'FontWeight', 'bold');
if n_pass > 0
    legend(ax3, 'Location', 'northeast', 'TextColor', txtCol, ...
        'Color', [0.08 0.10 0.14], 'EdgeColor', [0.30 0.32 0.38]);
end

%% Panel 4: Link margin vs elevation angle (scatter, colored by pass)
ax4 = subplot(2, 2, 4, 'Parent', fig);
set(ax4, 'Color', axCol, 'XColor', txtCol*0.8, 'YColor', txtCol*0.8, ...
    'GridColor', gridC, 'GridAlpha', 0.5);
hold(ax4, 'on');  grid(ax4, 'on');

for w = 1:n_pass
    col  = pass_cmap(min(w, size(pass_cmap,1)), :);
    mask = result.t_s >= passes(w).start_s & result.t_s <= passes(w).stop_s;
    if any(mask)
        scatter(ax4, result.el_deg(mask), result.LinkMargin_dB(mask), 12, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.7, 'DisplayName', sprintf('Pass %d', w));
    end
end

yline(ax4, 0, '--', 'Color', [0.95 0.35 0.35], 'LineWidth', 1.0, 'HandleVisibility', 'off');
xlabel(ax4, 'Elevation Angle (deg)', 'Color', txtCol, 'FontSize', 10);
ylabel(ax4, 'Link Margin (dB)',      'Color', txtCol, 'FontSize', 10);
title(ax4,  'Link Margin vs Elevation', 'Color', txtCol, 'FontSize', 10, 'FontWeight', 'bold');
if n_pass > 0
    legend(ax4, 'Location', 'southeast', 'TextColor', txtCol, ...
        'Color', [0.08 0.10 0.14], 'EdgeColor', [0.30 0.32 0.38]);
end

set(fig, 'Color', bgCol);
end
