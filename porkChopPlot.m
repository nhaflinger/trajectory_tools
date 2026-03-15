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

% Convert JD axis to calendar dates (approximate)
xt = linspace(min(departJD), max(departJD), 8);
xtLabels = datestr(xt - 1721058.5, 'yyyy-mm-dd');
set(gca, 'XTick', xt, 'XTickLabel', xtLabels);

title(sprintf('Pork Chop Plot: %s -> %s', departBody.name, arrivalBody.name));

end
