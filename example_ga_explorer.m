% example_ga_explorer.m
% Gravity-assist scenario explorer — Earth to Jupiter
%
% Evaluates six flyby architectures over a shared 2026-2032 launch window
% and compares their delta-V budgets.  Uses:
%
%   tisserandGraph    — design-space topology with iso-v∞ contours
%   resonantOrbits    — resonant return-orbit markers on the Tisserand graph
%   porkChopSequence  — per-leg launch-window overview for each sequence
%   flybySequence     — full patched-conic trajectory and ΔV budget
%
% Sequences evaluated:
%   1. Direct  E→J           baseline
%   2. E→V→J                 Venus gravity assist
%   3. E→M→J                 Mars gravity assist
%   4. E→E→J                 Earth resonance flyby
%   5. E→V→E→J               Venus + Earth flybys (VEEGA-style)
%   6. E→M→E→J  (MEGA)       Mars + Earth flybys (Europa Clipper style)
%                            M→E leg uses type2 (prograde long arc, >180°)
%
% Figures produced:
%   1  Tisserand graph (Venus/Earth/Mars/Jupiter) with resonant-orbit
%      markers and leg-orbit points for every computed sequence
%   2  Per-leg pork-chop for the lowest-DV single-flyby sequence
%   3  Per-leg pork-chop for the lowest-DV multi-flyby sequence
%   4  3D heliocentric trajectory for the lowest-DV sequence overall
%
% Runtime note: reduce nDepDates / nTofPoints in each searchOpts struct
% to trade resolution for speed.  Defaults give ~1-3 min total.

bodies = constants();
muSun  = bodies.Sun.mu;
AU     = bodies.Constants.AU;

outDir = fullfile(fileparts(mfilename('fullpath')), 'output_ga_explorer');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

jdStart = julianDate(2026,  1,  1);
jdEnd   = julianDate(2032, 12, 31);

% ---- Shared arrival / departure parameters ----
baseOpts = struct( ...
    'departureAltitude',     200,     ...   % km LEO parking orbit
    'departureInclination',  28.5,    ...   % deg KSC latitude
    'arrivalAltitude',       500,     ...   % km perijove
    'arrivalApogeeAltitude', 8100000, ...   % km highly elliptical capture (~Juno)
    'arrivalInclination',    0);

%% ---- Sequence definitions ------------------------------------------------

seqs = {};

% 1. Direct E→J (baseline) -------------------------------------------------
s             = struct();
s.name        = 'Direct  E->J';
s.bSeq        = {bodies.Earth, bodies.Jupiter};
s.tofRanges   = [600 1400];
s.seqOpts     = baseOpts;
s.searchOpts  = struct('nDepDates', 50, 'nTofPoints', 30);
seqs{end+1}   = s;

% 2. E→V→J -----------------------------------------------------------------
s             = struct();
s.name        = 'E->V->J';
s.bSeq        = {bodies.Earth, bodies.Venus, bodies.Jupiter};
s.tofRanges   = [100 250; 600 1200];
s.seqOpts     = baseOpts;
s.seqOpts.flybyAltitudes      = 300;
s.seqOpts.atmosphereAltitudes = 250;    % Venus atmosphere top ~250 km
s.searchOpts  = struct('nDepDates', 50, 'nTofPoints', 20, ...
    'flybyAltitudes', 300);
seqs{end+1}   = s;

% 3. E→M→J -----------------------------------------------------------------
s             = struct();
s.name        = 'E->M->J';
s.bSeq        = {bodies.Earth, bodies.Mars, bodies.Jupiter};
s.tofRanges   = [90 300; 600 1400];
s.seqOpts     = baseOpts;
s.seqOpts.flybyAltitudes      = 300;
s.seqOpts.atmosphereAltitudes = 50;     % Mars atmosphere top ~50 km
s.searchOpts  = struct('nDepDates', 50, 'nTofPoints', 20, ...
    'flybyAltitudes', 300);
seqs{end+1}   = s;

% 4. E→E→J  (Earth resonance flyby) ----------------------------------------
% Spacecraft departs Earth, completes a resonant heliocentric orbit
% (e.g. 2:3 ≈ 548 d, 1:2 ≈ 730 d), returns to Earth for a free flyby,
% then transfers to Jupiter.
s             = struct();
s.name        = 'E->E->J';
s.bSeq        = {bodies.Earth, bodies.Earth, bodies.Jupiter};
s.tofRanges   = [300 730; 600 1400];
s.seqOpts     = baseOpts;
s.seqOpts.flybyAltitudes      = 500;
s.seqOpts.atmosphereAltitudes = 100;    % Earth atmosphere top ~100 km
s.searchOpts  = struct('nDepDates', 50, 'nTofPoints', 20, ...
    'flybyAltitudes', 500);
seqs{end+1}   = s;

% 5. E→V→E→J  (Venus + Earth flybys, VEEGA-style) --------------------------
s             = struct();
s.name        = 'E->V->E->J';
s.bSeq        = {bodies.Earth, bodies.Venus, bodies.Earth, bodies.Jupiter};
s.tofRanges   = [100 250; 300 700; 600 1400];
s.seqOpts     = baseOpts;
s.seqOpts.flybyAltitudes      = [300 500];
s.seqOpts.atmosphereAltitudes = [250 100];
s.searchOpts  = struct('nDepDates', 35, 'nTofPoints', 12, ...
    'flybyAltitudes', [300 500]);
seqs{end+1}   = s;

% 6. E→M→E→J  (MEGA — Mars + Earth flybys, Europa Clipper style) -----------
% M→E leg is type2: CCW angle from Mars to Earth is >180° for the
% prograde long arc that sweeps around before the Earth flyby.
s             = struct();
s.name        = 'E->M->E->J (MEGA)';
s.bSeq        = {bodies.Earth, bodies.Mars, bodies.Earth, bodies.Jupiter};
s.tofRanges   = [90 200; 400 900; 700 1400];
s.seqOpts     = baseOpts;
s.seqOpts.flybyAltitudes      = [300 500];
s.seqOpts.atmosphereAltitudes = [50  100];
s.seqOpts.transferTypes       = {'type1', 'type2', 'type1'};
s.searchOpts  = struct('nDepDates', 35, 'nTofPoints', 12, ...
    'flybyAltitudes', [300 500]);
s.searchOpts.transferTypes    = {'type1', 'type2', 'type1'};
seqs{end+1}   = s;

nSeqs = numel(seqs);

%% ---- Grid search for each sequence ---------------------------------------

fprintf('Earth -> Jupiter  gravity-assist scenario search  (%s to %s)\n', ...
    datestr(jdStart - 1721058.5,'yyyy'), datestr(jdEnd - 1721058.5,'yyyy'));
fprintf('%s\n\n', repmat('=',1,62));

bestGrid    = cell(nSeqs, 1);
fullResults = cell(nSeqs, 1);

for i = 1:nSeqs
    s = seqs{i};
    fprintf('[%d/%d]  %-22s  searching ...', i, nSeqs, s.name);
    try
        best = findBestFlybyWindow(s.bSeq, jdStart, jdEnd, ...
            s.tofRanges, s.searchOpts);
        res  = flybySequence(s.bSeq, best.departureJD, best.tofDays, s.seqOpts);
        bestGrid{i}    = best;
        fullResults{i} = res;
        fprintf('  dep %s   DV %.3f km/s   TOF %.0f d\n', ...
            datestr(best.departureJD - 1721058.5,'yyyy-mmm-dd'), ...
            res.deltaV, res.tof);
    catch err
        fprintf('  FAILED: %s\n', err.message);
    end
end

%% ---- Comparison table ----------------------------------------------------

fprintf('\n%-24s | %8s | %8s | %8s | %9s | %12s\n', ...
    'Sequence', 'Dep DV', 'GA DV', 'Arr DV', 'Total DV', 'TOF');
fprintf('%-24s | %8s | %8s | %8s | %9s | %12s\n', ...
    '', '(km/s)', '(km/s)', '(km/s)', '(km/s)', '(days / yr)');
fprintf('%s\n', repmat('-',1,82));

for i = 1:nSeqs
    if isempty(fullResults{i}), continue; end
    r = fullResults{i};
    d = r.details;
    ga_str = '';
    if d.dvPoweredFlybys > 0.001
        ga_str = sprintf('%.3f', d.dvPoweredFlybys);
    else
        ga_str = 'free';
    end
    fprintf('%-24s | %8.3f | %8s | %8.3f | %9.3f | %6.0f / %.1f\n', ...
        seqs{i}.name, d.dvDeparture, ga_str, d.dvArrival, ...
        r.deltaV, r.tof, r.tof/365.25);
end
fprintf('%s\n', repmat('-',1,82));

% Best overall
bestTotalDV = Inf;
bestIdx     = 1;
for i = 1:nSeqs
    if ~isempty(fullResults{i}) && fullResults{i}.deltaV < bestTotalDV
        bestTotalDV = fullResults{i}.deltaV;
        bestIdx     = i;
    end
end
fprintf('\nLowest total DV: %s  (%.3f km/s)\n\n', ...
    seqs{bestIdx}.name, bestTotalDV);

%% ---- Resonant orbit tables -----------------------------------------------

resonantOrbits(bodies.Venus);
resonantOrbits(bodies.Earth);
resonantOrbits(bodies.Mars);

%% ---- Figure 1: Tisserand graph + resonant orbits + trajectory overlays ---

tBodies    = {bodies.Venus, bodies.Earth, bodies.Mars, bodies.Jupiter};
[fig1, ax1] = tisserandGraph(tBodies);

% Mark resonant return orbits for key flyby bodies (suppress table reprint)
resonantOrbits(bodies.Venus, 'ax', ax1, 'print', false);
resonantOrbits(bodies.Earth, 'ax', ax1, 'print', false);
resonantOrbits(bodies.Mars,  'ax', ax1, 'print', false);

% Overlay each sequence's transfer-ellipse (rp,ra) points.
% Each leg maps to one point: departure of that leg's Lambert arc.
% Gravity-assist moves the point along one body's v∞ contour.
seqColors = lines(nSeqs);
hold(ax1, 'on');
for i = 1:nSeqs
    if isempty(fullResults{i}), continue; end
    r     = fullResults{i};
    nLegs = numel(r.legs);
    rpAU  = nan(1, nLegs);
    raAU  = nan(1, nLegs);
    for k = 1:nLegs
        [rp, ra] = localRpRa(r.legs(k).r1_vec, r.legs(k).v1_transfer, muSun);
        rpAU(k)  = rp / AU;
        raAU(k)  = ra / AU;
    end
    valid = isfinite(rpAU) & isfinite(raAU);
    col   = seqColors(i,:);
    plot(ax1, rpAU(valid), raAU(valid), 'o--', ...
        'Color', col, 'MarkerSize', 7, 'MarkerFaceColor', col, ...
        'LineWidth', 1.3, 'DisplayName', seqs{i}.name);
end
% Refresh legend to include sequence overlays
legend(ax1, 'show', 'Location', 'northwest', 'FontSize', 7, ...
    'TextColor', [0.85 0.85 0.85], 'Color', [0.08 0.08 0.12], ...
    'EdgeColor', [0.4 0.4 0.4]);
saveas(fig1, fullfile(outDir, 'tisserand_graph.png'));

%% ---- Figure 2: Pork-chop for best single-flyby sequence ------------------

best1FB_dv  = Inf;
best1FB_idx = 2;   % default to E→V→J
for i = 1:nSeqs
    if numel(seqs{i}.bSeq) == 3 && ~isempty(fullResults{i})
        if fullResults{i}.deltaV < best1FB_dv
            best1FB_dv  = fullResults{i}.deltaV;
            best1FB_idx = i;
        end
    end
end

s      = seqs{best1FB_idx};
pcOpts = struct('markBest', bestGrid{best1FB_idx});
if isfield(s.seqOpts, 'transferTypes')
    pcOpts.transferTypes = s.seqOpts.transferTypes;
end
fig2 = porkChopSequence(s.bSeq, jdStart, jdEnd, s.tofRanges, pcOpts);
sgtitle(fig2, sprintf('Best 1-flyby: %s  (%.3f km/s total ΔV)', ...
    s.name, fullResults{best1FB_idx}.deltaV), ...
    'Color',[0.9 0.9 0.9], 'FontSize',11);
saveas(fig2, fullfile(outDir, 'pork_chop_best_single_flyby.png'));

%% ---- Figure 3: Pork-chop for best multi-flyby sequence -------------------

best2FB_dv  = Inf;
best2FB_idx = nSeqs;   % default to MEGA
for i = 1:nSeqs
    if numel(seqs{i}.bSeq) >= 4 && ~isempty(fullResults{i})
        if fullResults{i}.deltaV < best2FB_dv
            best2FB_dv  = fullResults{i}.deltaV;
            best2FB_idx = i;
        end
    end
end

s       = seqs{best2FB_idx};
pcOpts2 = struct('markBest', bestGrid{best2FB_idx});
if isfield(s.seqOpts, 'transferTypes')
    pcOpts2.transferTypes = s.seqOpts.transferTypes;
end
fig3 = porkChopSequence(s.bSeq, jdStart, jdEnd, s.tofRanges, pcOpts2);
sgtitle(fig3, sprintf('Best multi-flyby: %s  (%.3f km/s total ΔV)', ...
    s.name, fullResults{best2FB_idx}.deltaV), ...
    'Color',[0.9 0.9 0.9], 'FontSize',11);
saveas(fig3, fullfile(outDir, 'pork_chop_best_multi_flyby.png'));

%% ---- Figure 4: 3D trajectory for lowest-DV sequence ---------------------

fprintf('\nPlotting 3D trajectory for best sequence: %s\n', seqs{bestIdx}.name);
plotFlybySequence(fullResults{bestIdx}, seqs{bestIdx}.bSeq);
saveas(gcf, fullfile(outDir, '3d_trajectory_best_sequence.png'));

%% ---- Local helper --------------------------------------------------------

function [rp, ra] = localRpRa(r_vec, v_vec, mu)
%LOCALRPRA  Periapsis and apoapsis of the Keplerian orbit defined by (r,v).
r_vec = r_vec(:);  v_vec = v_vec(:);
R  = norm(r_vec);
a  = 1 / (2/R - dot(v_vec,v_vec)/mu);
if ~isfinite(a) || a <= 0
    rp = NaN;  ra = NaN;  return;
end
h_vec = cross(r_vec, v_vec);
e_vec = cross(v_vec, h_vec)/mu - r_vec/R;
e     = norm(e_vec);
if e >= 1.0
    rp = NaN;  ra = NaN;  return;
end
rp = a*(1 - e);
ra = a*(1 + e);
end
