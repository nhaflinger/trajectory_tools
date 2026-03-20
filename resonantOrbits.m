function result = resonantOrbits(body, varargin)
%RESONANTORBITS  Resonant return-orbit analysis for gravity-assist sequencing.
%
%   result = resonantOrbits(body)
%   result = resonantOrbits(body, options)
%   result = resonantOrbits(body, 'ax', axHandle)
%
%   For a p:q resonant orbit (spacecraft completes p revolutions while the
%   planet completes q), computes the required semi-major axis, orbital
%   period, minimum v∞ at a tangent encounter, and the apsis distances
%   assuming the flyby occurs at periapsis (outward orbits) or apoapsis
%   (inward orbits).
%
%   Physics:
%     a_sc  = a_p * (q/p)^(2/3)                 [Kepler's 3rd law]
%     v∞_min = v_p * |sqrt(2 - a_p/a_sc) - 1|   [tangent flyby, min energy]
%
%   Inputs:
%     body        - body struct from constants() — the flyby planet
%     options     (optional struct or name-value pairs)
%       .maxN      max integer in p and q numerator/denominator  [5]
%       .ax        axes handle — marks resonances on Tisserand graph
%       .print     print summary table                           [true]
%
%   Output:
%     result  - struct array, one element per resonance, with fields:
%       .label    string 'p:q'
%       .p, .q    numerator / denominator integers
%       .a_sc     spacecraft semi-major axis (km)
%       .T_days   spacecraft orbital period (days)
%       .v_inf    minimum v∞ at tangent encounter (km/s)
%       .r_p      periapsis of resonant orbit (km)
%       .r_a      apoapsis of resonant orbit (km)

% ---- parse variable arguments ----
opts   = struct('maxN', 5, 'print', true);
axHandle = [];
k = 1;
while k <= numel(varargin)
    arg = varargin{k};
    if ischar(arg) && strcmpi(arg,'ax') && k < numel(varargin)
        axHandle = varargin{k+1};  k = k+2;
    elseif ischar(arg) && strcmpi(arg,'maxN') && k < numel(varargin)
        opts.maxN = varargin{k+1};  k = k+2;
    elseif ischar(arg) && strcmpi(arg,'print') && k < numel(varargin)
        opts.print = varargin{k+1};  k = k+2;
    elseif isstruct(arg)
        for fi = fieldnames(arg)', opts.(fi{1}) = arg.(fi{1}); end
        k = k+1;
    else
        k = k+1;
    end
end

consts = constants();
muSun  = consts.Sun.mu;
AU     = consts.Constants.AU;

a_p   = body.a;
T_p   = 2*pi * sqrt(a_p^3 / muSun);   % s
T_p_d = T_p / 86400;                   % days
v_p   = sqrt(muSun / a_p);            % km/s

maxN = opts.maxN;

% ---- enumerate reduced-fraction resonances p:q ----
pairs = zeros(0,2);
for pp = 1:maxN
    for qq = 1:maxN
        if pp ~= qq && gcd(pp,qq) == 1
            pairs(end+1,:) = [pp, qq]; %#ok<AGROW>
        end
    end
end
% Sort by period ratio q/p (ascending: inward first, then outward)
[~,idx] = sort(pairs(:,2) ./ pairs(:,1));
pairs   = pairs(idx,:);

% ---- compute resonance parameters ----
result = struct('label',{},'p',{},'q',{},'a_sc',{},'T_days',{}, ...
    'v_inf',{},'r_p',{},'r_a',{});

for i = 1:size(pairs,1)
    pp = pairs(i,1);  qq = pairs(i,2);
    a_sc = a_p * (qq/pp)^(2/3);

    % Minimum v∞ (tangent encounter — spacecraft and planet velocities parallel)
    v_inf = v_p * abs(sqrt(max(0, 2 - a_p/a_sc)) - 1);

    % Apsis distances assuming flyby at periapsis (outward) or apoapsis (inward)
    if a_sc >= a_p
        r_p = a_p;
        r_a = 2*a_sc - a_p;
    else
        r_a = a_p;
        r_p = 2*a_sc - a_p;
        if r_p <= 0
            continue;   % orbit doesn't extend to planet — not a valid resonance
        end
    end

    % Skip if periapsis is inside the body
    if r_p < body.radius * 1.05
        continue;
    end

    entry.label  = sprintf('%d:%d', pp, qq);
    entry.p      = pp;
    entry.q      = qq;
    entry.a_sc   = a_sc;
    entry.T_days = T_p_d * (qq/pp);
    entry.v_inf  = v_inf;
    entry.r_p    = r_p;
    entry.r_a    = r_a;

    if isempty(result)
        result = entry;
    else
        result(end+1) = entry; %#ok<AGROW>
    end
end

% ---- print table ----
if opts.print
    fprintf('\nResonant Return Orbits at %s  (a_p = %.4f AU, T_p = %.1f days)\n', ...
        body.name, a_p/AU, T_p_d);
    fprintf('%s\n', repmat('-',1,73));
    fprintf('  Ratio  |  T_sc (d) | a_sc (AU) | v_inf_min (km/s) |  r_p (AU) |  r_a (AU)\n');
    fprintf('%s\n', repmat('-',1,73));
    for i = 1:numel(result)
        r = result(i);
        dir = '';
        if r.a_sc > a_p, dir = 'out'; else, dir = ' in'; end
        fprintf('  %5s  | %9.2f | %9.5f |        %8.4f   | %9.5f | %9.5f  %s\n', ...
            r.label, r.T_days, r.a_sc/AU, r.v_inf, r.r_p/AU, r.r_a/AU, dir);
    end
    fprintf('%s\n\n', repmat('-',1,73));
    fprintf('  v_inf_min = v_p * |sqrt(2 - a_p/a_sc) - 1|  (tangent encounter)\n');
    fprintf('  r_p/r_a computed with flyby at periapsis (outward) or apoapsis (inward)\n\n');
end

% ---- optional Tisserand graph overlay ----
if ~isempty(axHandle)
    hold(axHandle,'on');
    col = bodyColor(body.name);
    for i = 1:numel(result)
        r = result(i);
        rp_AU = r.r_p / AU;
        ra_AU = r.r_a / AU;
        plot(axHandle, rp_AU, ra_AU, 'x', 'MarkerSize',9, ...
            'Color',col, 'LineWidth',1.8, 'HandleVisibility','off');
        text(axHandle, rp_AU*1.02, ra_AU*1.04, ...
            sprintf(' %s\n %.1f km/s', r.label, r.v_inf), ...
            'Color',col,'FontSize',6.5,'HandleVisibility','off');
    end
end
end

% =========================================================================
function col = bodyColor(name)
switch lower(name)
    case 'mercury',  col = [0.60 0.58 0.55];
    case 'venus',    col = [0.85 0.75 0.35];
    case 'earth',    col = [0.20 0.45 0.75];
    case 'mars',     col = [0.72 0.28 0.18];
    case 'jupiter',  col = [0.75 0.60 0.45];
    case 'saturn',   col = [0.85 0.80 0.55];
    case 'uranus',   col = [0.55 0.85 0.85];
    case 'neptune',  col = [0.25 0.40 0.85];
    otherwise,       col = [0.70 0.70 0.70];
end
end
