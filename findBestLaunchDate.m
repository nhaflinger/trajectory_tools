function best = findBestLaunchDate(departBody, arrivalBody, jdStart, jdEnd, tofDays, options)
%FINDBESTLAUNCHDATE Find the best departure JD (min delta-v) over a window
%   best = findBestLaunchDate(departBody, arrivalBody, jdStart, jdEnd, tofDays, options)
%   Searches departure dates between jdStart and jdEnd (inclusive) and
%   time-of-flight values in tofDays (days) to find the minimum delta-v
%   transfer using a Lambert solver on circular orbits.
%
%   Output struct has fields:
%     .departureJD
%     .tofDays
%     .deltaV
%     .dvGrid (matrix)
%   and the input grids for reference.

if nargin < 6 || isempty(options)
    options = struct();
end
if ~isfield(options, 'epoch0') || isempty(options.epoch0)
    options.epoch0 = 2451545.0;
end
if ~isfield(options, 'mu') || isempty(options.mu)
    options.mu = constants().Sun.mu;
end

% Use 5-day steps if not provided
if numel(jdStart) == 1 && numel(jdEnd) == 1
    departJD = (jdStart:5:jdEnd)';
else
    departJD = jdStart;
end

tofDays = tofDays(:);

nd = numel(departJD);
tn = numel(tofDays);

dvGrid = nan(tn, nd);

for i = 1:nd
    jd = departJD(i);
    [r1, v1_circ] = orbitalStateCircular(departBody, jd, options.epoch0);
    for j = 1:tn
        tof = tofDays(j) * 86400;
        [r2, v2_circ] = orbitalStateCircular(arrivalBody, jd + tofDays(j), options.epoch0);
        try
            [v1t, v2t] = lambertSolver(r1, r2, tof, options.mu);
            dvGrid(j,i) = norm(v1t - v1_circ) + norm(v2_circ - v2t);
        catch
            dvGrid(j,i) = NaN;
        end
    end
end

[minDv, idx] = min(dvGrid(:));
[row, col] = ind2sub(size(dvGrid), idx);

best.departureJD = departJD(col);
best.tofDays = tofDays(row);
best.deltaV = minDv;
best.dvGrid = dvGrid;
best.departJD = departJD;
best.tofDaysGrid = tofDays;

end
