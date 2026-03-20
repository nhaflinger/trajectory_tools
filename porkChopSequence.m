function fig = porkChopSequence(bodies, jdStart, jdEnd, tofRanges, options)
%PORKCHOPSEQUENCE  Per-leg pork-chop plots for a multi-body flyby sequence.
%
%   fig = porkChopSequence(bodies, jdStart, jdEnd, tofRanges)
%   fig = porkChopSequence(bodies, jdStart, jdEnd, tofRanges, options)
%
%   Creates one subplot per leg showing hyperbolic excess speed (v∞) at
%   arrival as a function of departure date and time of flight.  Each leg
%   is solved independently via the Lambert problem.  The departure window
%   for leg i is automatically offset from jdStart/jdEnd by the accumulated
%   minimum TOF of all preceding legs.
%
%   Use the plots together to identify compatible windows: a free gravity
%   assist requires the v∞ at arrival of leg i ≈ v∞ at departure of leg i+1.
%
%   Inputs:
%     bodies      - 1×N cell array of body structs from constants()
%     jdStart     - Julian Date opening the leg-1 departure window
%     jdEnd       - Julian Date closing the leg-1 departure window
%     tofRanges   - (N-1)×2 matrix; tofRanges(i,:) = [min_days, max_days]
%     options     (optional struct)
%       .nDep          departure-date grid points per leg    [60]
%       .nTof          TOF grid points per leg               [60]
%       .transferTypes 1×(N-1) cell of 'type1'|'type2'       ['type1' each]
%       .markBest      struct with .departureJD + .tofDays to mark as ★
%       .cLimPct       colorscale percentile clamp [lo hi]   [5 85]
%
%   Output:
%     fig  - handle to the figure

if nargin < 5, options = struct(); end
if ~isfield(options,'nDep'),    options.nDep    = 60; end
if ~isfield(options,'nTof'),    options.nTof    = 60; end
if ~isfield(options,'cLimPct'), options.cLimPct = [5 85]; end

N     = numel(bodies);
nLegs = N - 1;

if ~isfield(options,'transferTypes')
    options.transferTypes = repmat({'type1'}, 1, nLegs);
end
if ischar(options.transferTypes)
    options.transferTypes = repmat({options.transferTypes}, 1, nLegs);
end

muSun = constants().Sun.mu;

% ---- figure layout ----
ncols = min(nLegs, 2);
nrows = ceil(nLegs / ncols);
fig = figure('Name', sprintf('Pork-Chop Sequence: %s', seqName(bodies)), ...
    'NumberTitle','off', 'Color',[0.06 0.06 0.10]);

for iLeg = 1:nLegs
    % Departure JD window for this leg (offset by min TOF of prior legs)
    shift_lo = sum(tofRanges(1:iLeg-1, 1));
    shift_hi = sum(tofRanges(1:iLeg-1, 2));
    jd_lo = jdStart + shift_lo;
    jd_hi = jdEnd   + shift_hi;

    depJD_vec = linspace(jd_lo, jd_hi, options.nDep);
    tofVec    = linspace(tofRanges(iLeg,1), tofRanges(iLeg,2), options.nTof);
    isLW      = strcmpi(options.transferTypes{iLeg}, 'type2');

    bDep = bodies{iLeg};
    bArr = bodies{iLeg+1};

    % ---- compute v∞ at arrival on the grid ----
    vInfGrid = nan(options.nTof, options.nDep);
    for id = 1:options.nDep
        jd_dep = depJD_vec(id);
        [r1, ~] = orbitalState(bDep, jd_dep);
        for it = 1:options.nTof
            jd_arr = jd_dep + tofVec(it);
            [r2, v2_body] = orbitalState(bArr, jd_arr);
            try
                [~, v2t] = lambertSolver(r1, r2, tofVec(it)*86400, muSun, isLW);
                vInfGrid(it, id) = norm(v2t - v2_body);
            catch
            end
        end
    end

    % ---- subplot ----
    figure(fig);
    ax = subplot(nrows, ncols, iLeg);
    set(ax, 'Color',[0.08 0.08 0.12], ...
        'XColor',[0.7 0.7 0.7], 'YColor',[0.7 0.7 0.7]);

    % Colorscale — clamp to percentile range so minimum-energy island is visible
    vFlat = vInfGrid(isfinite(vInfGrid));
    if numel(vFlat) > 2
        clo = localPrctile(vFlat, options.cLimPct(1));
        chi = localPrctile(vFlat, options.cLimPct(2));
    else
        clo = 0;  chi = 20;
    end

    [~,~] = contourf(ax, depJD_vec, tofVec, vInfGrid, 20, 'LineColor','none');
    colormap(ax, parula);
    caxis(ax, [clo chi]);
    cb = colorbar(ax);
    cb.Label.String = 'v_{\infty} arrive (km/s)';
    cb.Color = [0.7 0.7 0.7];

    hold(ax,'on');

    % ---- mark best point ----
    if isfield(options,'markBest')
        mb = options.markBest;
        if isfield(mb,'departureJD') && isfield(mb,'tofDays') && ...
                numel(mb.tofDays) >= iLeg
            jd_best = mb.departureJD;
            for k = 1:iLeg-1
                jd_best = jd_best + mb.tofDays(k);
            end
            plot(ax, jd_best, mb.tofDays(iLeg), 'p', 'MarkerSize',13, ...
                'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'LineWidth',1.0);
        end
    end

    % ---- labels / ticks ----
    typStr = options.transferTypes{iLeg};
    title(ax, sprintf('Leg %d: %s \\rightarrow %s  [%s]', ...
        iLeg, bDep.name, bArr.name, typStr), 'Color',[0.9 0.9 0.9]);
    ylabel(ax, 'TOF (days)', 'Color',[0.8 0.8 0.8]);

    setCalendarTicks(ax, depJD_vec(1), depJD_vec(end));
end

sgtitle(fig, sprintf('Pork-Chop Sequence: %s', seqName(bodies)), ...
    'Color',[0.9 0.9 0.9], 'FontSize',11);
end

% =========================================================================
%  Local helpers
% =========================================================================

function setCalendarTicks(ax, jd_lo, jd_hi)
dn_min = jd_lo - 1721058.5;
dn_max = jd_hi - 1721058.5;
dv_lo  = datevec(dn_min);
dv_hi  = datevec(dn_max);

% Major ticks: Jan 1 of each year in range
yr_dns = arrayfun(@(y) datenum(y,1,1), dv_lo(1):dv_hi(1)+1);
yr_dns = yr_dns(yr_dns >= dn_min & yr_dns <= dn_max);

% Minor ticks: 1st of each non-January month
y = dv_lo(1);  m = dv_lo(2)+1;
if m > 12, m=1; y=y+1; end
mo_dns = [];
while datenum(y,m,1) <= dn_max
    if m ~= 1, mo_dns(end+1) = datenum(y,m,1); end %#ok<AGROW>
    m = m+1;
    if m > 12, m=1; y=y+1; end
end

if numel(yr_dns) >= 2
    set(ax, 'XTick', yr_dns+1721058.5, 'XTickLabel', datestr(yr_dns(:),'yyyy'));
    if ~isempty(mo_dns)
        ax.XAxis.MinorTickValues = mo_dns + 1721058.5;
        ax.XMinorTick = 'on';
    end
else
    all_dns = sort([yr_dns, mo_dns]);
    if ~isempty(all_dns)
        set(ax,'XTick', all_dns+1721058.5, 'XTickLabel', datestr(all_dns(:),'mmm yyyy'));
        xtickangle(ax, 45);
    end
end
ax.XLabel.String = 'Departure Date';
ax.XLabel.Color  = [0.8 0.8 0.8];
end

% -------------------------------------------------------------------------
function p = localPrctile(x, pct)
% Percentile without Statistics Toolbox
x   = sort(x(isfinite(x)));
n   = numel(x);
if n == 0, p = 0; return; end
idx = max(1, min(n, round(pct/100 * n)));
p   = x(idx);
end

% -------------------------------------------------------------------------
function s = seqName(bodies)
names = cellfun(@(b) b.name, bodies, 'UniformOutput',false);
s = strjoin(names,' -> ');
end
