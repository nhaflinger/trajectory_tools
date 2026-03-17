function porkChopPlot(departBody, arrivalBody, departJD, tofDays, options)
%PORKCHOPPLOT Generate a basic pork-chop plot (departure date vs time-of-flight)
%   porkChopPlot(departBody, arrivalBody, departJD, tofDays, options)
%
%   departJD is a vector of Julian Dates for departure.
%   tofDays is a vector of time-of-flight values in days.
%   options can contain:
%     .epoch0 - reference Julian date for orbital phasing (default J2000)
%     .mu     - central body gravitational parameter (default Sun mu)
%     .fixPlane - whether to assume coplanar orbits (default true)
%
%   The function computes an approximate delta-V for each (depart, tof)
%   combination using a Lambert solution on circular orbits.

% Backwards compatibility for MATLAB versions without arguments blocks
if nargin < 5 || isempty(options)
    options = struct();
end
if ~isfield(options, 'epoch0') || isempty(options.epoch0)
    options.epoch0 = 2451545.0;
end
if ~isfield(options, 'mu') || isempty(options.mu)
    options.mu = constants().Sun.mu;
end
if ~isfield(options, 'fixPlane') || isempty(options.fixPlane)
    options.fixPlane = true;
end

nd = numel(departJD);
tn = numel(tofDays);

dvGrid = nan(tn, nd);

for i = 1:nd
    jd = departJD(i);
    [r1, v1_circ] = orbitalStateCircular(departBody, jd, options.epoch0);
    for j = 1:tn
        tof = tofDays(j) * 86400;
        % propagate arrival body to arrival epoch
        [r2, v2_circ] = orbitalStateCircular(arrivalBody, jd + tofDays(j), options.epoch0);
        try
            [v1t, v2t] = lambertSolver(r1, r2, tof, options.mu);
            dvGrid(j,i) = norm(v1t - v1_circ) + norm(v2_circ - v2t);
        catch
            dvGrid(j,i) = NaN;
        end
    end
end

% Plot
figure('Name', sprintf('Pork Chop: %s -> %s', departBody.name, arrivalBody.name), 'NumberTitle', 'off');
[TT, DD] = meshgrid(departJD, tofDays);
contourf(TT, DD, dvGrid, 20, 'LineColor','none');
colormap(parula);
cb = colorbar;
cb.Label.String = 'ΔV (km/s)';
xlabel('Departure Date');
ylabel('Time of Flight (days)');

% Dual-level calendar ticks: year boundaries (major) + month starts (minor).
% For ranges < 1 year with no Jan 1 present, falls back to monthly major ticks.
dn_min = min(departJD) - 1721058.5;   % JD -> MATLAB datenum
dn_max = max(departJD) - 1721058.5;
ax     = gca;

% --- Major ticks: Jan 1 of each year within range ---
dv_lo   = datevec(dn_min);
dv_hi   = datevec(dn_max);
yr_dns  = arrayfun(@(y) datenum(y, 1, 1), dv_lo(1) : dv_hi(1)+1);
yr_dns  = yr_dns(yr_dns >= dn_min & yr_dns <= dn_max);

% --- Minor ticks: 1st of each non-January month within range ---
y = dv_lo(1);  m = dv_lo(2) + 1;
if m > 12,  m = 1;  y = y + 1;  end
mo_dns = [];
while datenum(y, m, 1) <= dn_max
    if m ~= 1
        mo_dns(end+1) = datenum(y, m, 1);  %#ok<AGROW>
    end
    m = m + 1;
    if m > 12,  m = 1;  y = y + 1;  end
end

if numel(yr_dns) >= 2
    % Multi-year span: year labels on major ticks, unlabelled minor month ticks
    set(ax, 'XTick', yr_dns + 1721058.5, 'XTickLabel', datestr(yr_dns(:), 'yyyy'));
    if ~isempty(mo_dns)
        ax.XAxis.MinorTickValues = mo_dns + 1721058.5;
        ax.XMinorTick = 'on';
    end
else
    % Sub-year span: monthly major ticks with month+year labels
    all_dns = sort([yr_dns, mo_dns]);
    set(ax, 'XTick', all_dns + 1721058.5, 'XTickLabel', datestr(all_dns(:), 'mmm yyyy'));
    xtickangle(45);
end

title(sprintf('Pork Chop Plot: %s -> %s', departBody.name, arrivalBody.name));

end
