function best = findBestFlybyWindow(bodies, jdStart, jdEnd, tofRanges, options)
%FINDBESTFLYBYWWINDOW  Grid search for minimum-ΔV gravity-assist trajectory.
%
%   best = findBestFlybyWindow(bodies, jdStart, jdEnd, tofRanges)
%   best = findBestFlybyWindow(bodies, jdStart, jdEnd, tofRanges, options)
%
%   Searches over departure dates and per-leg TOFs to minimise the total
%   heliocentric ΔV proxy: v∞_departure + Σ(powered-flyby ΔV) + v∞_arrival.
%   Departure and arrival burns are body-centric and computed in post-
%   processing by flybySequence; the grid search uses the heliocentric
%   v∞ magnitudes as a fast objective function.
%
%   Inputs:
%     bodies      - 1×N cell array of body structs from constants()
%                   (first = departure, last = arrival, middle = flyby)
%     jdStart     - earliest departure Julian Date
%     jdEnd       - latest departure Julian Date
%     tofRanges   - (N−1)×2 matrix: [min_days, max_days] per leg
%                   OR (N−1)×M matrix of explicit TOF values (rows = legs)
%     options
%       .nDepDates      number of departure dates to sample [50]
%       .nTofPoints     number of TOF samples per leg [25]
%       .flybyAltitudes 1×(N−2) powered-flyby periapsis altitudes (km) [300]
%       .transferTypes  per-leg 'type1'|'type2' cell or scalar ['type1']
%
%   Output (best struct):
%     .departureJD  optimal departure Julian Date
%     .tofDays      1×(N−1) optimal TOFs per leg (days)
%     .deltaVProxy  heliocentric ΔV proxy at the best point (km/s)
%     .departJD     full departure-date grid used in the search
%     .tofGrids     cell array of TOF grids used per leg

if nargin < 5, options = struct(); end
if ~isfield(options, 'nDepDates'),  options.nDepDates  = 50; end
if ~isfield(options, 'nTofPoints'), options.nTofPoints = 25; end

N     = numel(bodies);
nLegs = N - 1;
nFB   = N - 2;

if ~isfield(options, 'flybyAltitudes')
    options.flybyAltitudes = 300 * ones(1, max(nFB, 1));
end
if ~isfield(options, 'transferTypes')
    options.transferTypes = repmat({'type1'}, 1, nLegs);
end
if ischar(options.transferTypes)
    options.transferTypes = repmat({options.transferTypes}, 1, nLegs);
end

muSun = constants().Sun.mu;

% --- Build TOF grids per leg ---
tofGrids = cell(nLegs, 1);
for i = 1:nLegs
    if size(tofRanges, 2) == 2
        tofGrids{i} = linspace(tofRanges(i,1), tofRanges(i,2), options.nTofPoints);
    else
        tofGrids{i} = tofRanges(i, :);
    end
end
depJDs    = linspace(jdStart, jdEnd, options.nDepDates);
tofSizes  = cellfun(@numel, tofGrids);
nTofTotal = prod(tofSizes);

% Precompute stride for multi-index decomposition
stride = ones(1, nLegs);
for k = 2:nLegs
    stride(k) = stride(k-1) * tofSizes(k-1);
end

% --- Grid search ---
bestDV    = Inf;
bestDepJD = NaN;
bestTofs  = NaN(1, nLegs);

for iDep = 1:numel(depJDs)
    jd0 = depJDs(iDep);

    for iTof = 1:nTofTotal
        % Decompose flat index -> per-leg TOF values
        tofs = zeros(1, nLegs);
        temp = iTof - 1;
        for k = 1:nLegs
            idx_k   = mod(floor(temp / stride(k)), tofSizes(k)) + 1;
            tofs(k) = tofGrids{k}(idx_k);
        end

        % Solve all legs; skip if any Lambert fails
        solved    = true;
        v1t_legs  = cell(nLegs, 1);
        v2t_legs  = cell(nLegs, 1);
        v1b_legs  = cell(nLegs, 1);
        v2b_legs  = cell(nLegs, 1);
        jdAccum   = jd0;

        for i = 1:nLegs
            jd_dep_i = jdAccum;
            jd_arr_i = jd_dep_i + tofs(i);
            tof_i    = tofs(i) * 86400;
            jdAccum  = jd_arr_i;

            [r1, v1b] = orbitalState(bodies{i},   jd_dep_i);
            [r2, v2b] = orbitalState(bodies{i+1}, jd_arr_i);

            try
                isLW = strcmpi(options.transferTypes{i}, 'type2');
                [v1t, v2t] = lambertSolver(r1, r2, tof_i, muSun, isLW);
            catch
                solved = false;
                break;
            end

            v1t_legs{i} = v1t;  v2t_legs{i} = v2t;
            v1b_legs{i} = v1b;  v2b_legs{i} = v2b;
        end

        if ~solved, continue; end

        % Heliocentric ΔV proxy
        v_inf_dep = norm(v1t_legs{1}     - v1b_legs{1});
        v_inf_arr = norm(v2b_legs{nLegs} - v2t_legs{nLegs});

        dv_powered = 0;
        for j = 1:nFB
            v_in  = norm(v2t_legs{j}   - v2b_legs{j});
            v_out = norm(v1t_legs{j+1} - v1b_legs{j+1});
            if abs(v_in - v_out) > 1e-4
                r_p     = bodies{j+1}.radius + options.flybyAltitudes(j);
                vp_in   = sqrt(v_in^2  + 2*bodies{j+1}.mu/r_p);
                vp_out  = sqrt(v_out^2 + 2*bodies{j+1}.mu/r_p);
                dv_powered = dv_powered + abs(vp_out - vp_in);
            end
        end

        dv_proxy = v_inf_dep + dv_powered + v_inf_arr;

        if dv_proxy < bestDV
            bestDV    = dv_proxy;
            bestDepJD = jd0;
            bestTofs  = tofs;
        end
    end
end

best.departureJD = bestDepJD;
best.tofDays     = bestTofs;
best.deltaVProxy = bestDV;
best.departJD    = depJDs;
best.tofGrids    = tofGrids;
end
