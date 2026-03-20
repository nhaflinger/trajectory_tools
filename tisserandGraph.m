function [fig, ax] = tisserandGraph(bodies, options)
%TISSERANDGRAPH  Tisserand parameter graph for gravity-assist mission design.
%
%   fig = tisserandGraph(bodies)
%   fig = tisserandGraph(bodies, options)
%
%   Plots iso-v∞ contours for each body on the (periapsis distance,
%   apoapsis distance) plane.  A free gravity assist moves the spacecraft
%   along one body's contour (same v∞ magnitude, deflected direction).
%   A deep-space manoeuvre moves it between contours.
%
%   Inputs:
%     bodies  - 1×N cell array of body structs from constants()
%     options (optional struct)
%       .vInfValues  v∞ contour levels in km/s        [0.5 1 2 3 5 10 20]
%       .nGrid       grid resolution                   [400]
%       .rpLim       periapsis range [min max] AU      [auto]
%       .raLim       apoapsis range  [min max] AU      [auto]
%       .trajectory  flybySequence result to overlay   (optional)
%       .showLabels  label contour values              [true]
%
%   Output:
%     fig  - handle to the figure

if nargin < 2, options = struct(); end
if ~isfield(options,'vInfValues'), options.vInfValues = [0.5 1 2 3 5 10 20]; end
if ~isfield(options,'nGrid'),      options.nGrid      = 400; end
if ~isfield(options,'showLabels'), options.showLabels = true; end

consts = constants();
muSun  = consts.Sun.mu;
AU     = consts.Constants.AU;

N    = numel(bodies);
aAU  = cellfun(@(b) b.a/AU, bodies);

if ~isfield(options,'rpLim')
    options.rpLim = [max(0.05, min(aAU)*0.5),  max(aAU)*1.08];
end
if ~isfield(options,'raLim')
    options.raLim = [min(aAU)*0.7,  max(aAU)*2.5];
end

nG    = options.nGrid;
rpVec = linspace(options.rpLim(1), options.rpLim(2), nG);
raVec = linspace(options.raLim(1), options.raLim(2), nG);
[RP_AU, RA_AU] = meshgrid(rpVec, raVec);
RP_km = RP_AU * AU;
RA_km = RA_AU * AU;

% ---- figure / axes ----
fig = figure('Name','Tisserand Graph','NumberTitle','off', ...
    'Color',[0.06 0.06 0.10]);
ax = axes('Parent',fig,'Color',[0.08 0.08 0.12], ...
    'XColor',[0.7 0.7 0.7],'YColor',[0.7 0.7 0.7]);
hold(ax,'on');  grid(ax,'on');
ax.GridColor = [0.25 0.25 0.25];
ax.GridAlpha = 0.5;

xlabel(ax,'Periapsis Distance (AU)','Color',[0.8 0.8 0.8]);
ylabel(ax,'Apoapsis Distance (AU)', 'Color',[0.8 0.8 0.8]);
title(ax,'Tisserand Graph — Gravity-Assist Design Space', ...
    'Color',[0.9 0.9 0.9]);

% Diagonal r_a = r_p (circular-orbit locus)
r_d = linspace(options.rpLim(1), options.rpLim(2), 300);
plot(ax, r_d, r_d, '--', 'Color',[0.5 0.5 0.6 0.5], 'LineWidth',0.8);
lbl_idx = round(numel(r_d)*0.72);
text(ax, r_d(lbl_idx)*1.01, r_d(lbl_idx)*1.06, 'circular', ...
    'Color',[0.5 0.5 0.6],'FontSize',7,'Rotation',35, ...
    'HorizontalAlignment','center');

% ---- v∞ contours per body ----
for i = 1:N
    b    = bodies{i};
    a_p  = b.a;               % km
    v_p  = sqrt(muSun / a_p); % km/s
    a_pAU = a_p / AU;

    % Tisserand parameter: T = a_p/a + 2*sqrt((a/a_p)*(1-e²))
    % 1-e² = 4*rp*ra / (rp+ra)²
    a_sc         = (RP_km + RA_km) / 2;
    one_minus_e2 = max(0, 4*RP_km.*RA_km ./ (RP_km + RA_km).^2);
    T            = a_p./a_sc + 2*sqrt((a_sc/a_p).*one_minus_e2);
    vInf         = real(sqrt(max(0, v_p^2*(3 - T))));

    % Valid region: orbit crosses planet distance AND ra ≥ rp
    valid = (RP_AU <= a_pAU) & (RA_AU >= a_pAU) & (RA_AU >= RP_AU);
    vInf(~valid) = NaN;

    col = bodyColor(b.name);

    [C,h] = contour(ax, rpVec, raVec, vInf, options.vInfValues, ...
        'Color',col,'LineWidth',1.1);
    h.HandleVisibility = 'off';   % managed in explicit legend below
    if options.showLabels && ~isempty(C)
        clabel(C, h, options.vInfValues, 'FontSize',6.5, 'Color',col, ...
            'LabelSpacing',300);
    end

    % Planet marker (circular orbit: rp = ra = a_planet)
    plot(ax, a_pAU, a_pAU, 'o', 'MarkerSize',9, ...
        'MarkerFaceColor',col, 'MarkerEdgeColor','w', 'LineWidth',0.8, ...
        'HandleVisibility','off');
    text(ax, a_pAU*1.03, a_pAU*1.07, ['  ' b.name], ...
        'Color',col,'FontSize',9,'FontWeight','bold');
end

% ---- optional trajectory overlay ----
if isfield(options,'trajectory')
    overlayTrajectory(ax, options.trajectory, muSun, AU);
end

xlim(ax, options.rpLim);
ylim(ax, options.raLim);

% ---- legend ----
% Collect explicit handles so every symbol type is explained.
legH   = {};
legTxt = {};

% Per-body: filled dot on a line represents planet position + v∞ contours
for i = 1:N
    col = bodyColor(bodies{i}.name);
    h = plot(ax, NaN, NaN, 'o-', 'Color',col, ...
        'MarkerFaceColor',col, 'MarkerEdgeColor','w', ...
        'MarkerSize',7, 'LineWidth',1.2);
    legH{end+1}   = h;                                          %#ok<AGROW>
    legTxt{end+1} = sprintf('%s  — planet (●) & v∞ contours', bodies{i}.name); %#ok<AGROW>
end

% Circular orbit locus (diagonal r_p = r_a)
h_circ = plot(ax, NaN, NaN, '--', 'Color',[0.5 0.5 0.6], 'LineWidth',0.8);
legH{end+1}   = h_circ;
legTxt{end+1} = 'Circular orbit  (r_p = r_a)';

% Resonant return orbits — × markers added by resonantOrbits(body, ''ax'', ax)
h_res = plot(ax, NaN, NaN, 'x', 'Color',[0.60 0.60 0.65], ...
    'MarkerSize',9, 'LineWidth',1.8);
legH{end+1}   = h_res;
legTxt{end+1} = '× p:q resonant return orbit  (label: ratio, min v∞)';

% Trajectory overlay entries (if a flybySequence result was supplied)
if isfield(options,'trajectory')
    h_path = plot(ax, NaN, NaN, 'w--', 'LineWidth',1.2);
    h_dot  = plot(ax, NaN, NaN, 'o', 'MarkerFaceColor',[0.55 0.55 0.55], ...
        'MarkerEdgeColor','w', 'MarkerSize',8, 'LineWidth',0.8);
    legH{end+1}   = h_path;
    legTxt{end+1} = 'Trajectory sequence path';
    legH{end+1}   = h_dot;
    legTxt{end+1} = 'L#  Transfer-leg orbit (r_p, r_a)';
end

legend(ax, [legH{:}], legTxt, 'Location','northwest', ...
    'TextColor',[0.85 0.85 0.85], 'Color',[0.10 0.10 0.15], ...
    'EdgeColor',[0.4 0.4 0.4], 'FontSize',8);
end

% =========================================================================
%  Local helpers
% =========================================================================

function overlayTrajectory(ax, result, muSun, AU)
%OVERLAYTRAJECTORY  Plot transfer-ellipse (rp, ra) points on Tisserand graph.
%
%   Each leg of the sequence maps to one point in the Tisserand plane.
%   A free gravity assist at body j moves the spacecraft from leg j's
%   arrival point to leg j+1's departure point along body j's v∞ contour.
nLegs = numel(result.legs);
legColors = [0.20 0.60 1.00;
             0.20 0.80 0.35;
             0.95 0.60 0.10;
             0.80 0.20 0.20;
             0.70 0.30 0.90];
if nLegs > size(legColors,1)
    legColors = repmat(legColors, ceil(nLegs/size(legColors,1)), 1);
end

rpArr = nan(1, nLegs);
raArr = nan(1, nLegs);
for i = 1:nLegs
    [rp, ra] = stateToRpRa(result.legs(i).r1_vec, result.legs(i).v1_transfer, muSun);
    rpArr(i) = rp / AU;
    raArr(i) = ra / AU;
end

% Dashed connecting line
plot(ax, rpArr, raArr, 'w--', 'LineWidth',1.2, 'HandleVisibility','off');

% Per-leg point
for i = 1:nLegs
    col = legColors(i,:);
    if isfinite(rpArr(i)) && isfinite(raArr(i))
        plot(ax, rpArr(i), raArr(i), 'o', 'MarkerSize',9, ...
            'MarkerFaceColor',col, 'MarkerEdgeColor','w', 'LineWidth',0.8, ...
            'HandleVisibility','off');
        text(ax, rpArr(i), raArr(i)*1.06, sprintf(' L%d',i), ...
            'Color',[0.9 0.9 0.9],'FontSize',7,'HandleVisibility','off');
    end
end
end

% -------------------------------------------------------------------------
function [rp, ra] = stateToRpRa(r_vec, v_vec, mu)
r_vec = r_vec(:);  v_vec = v_vec(:);
R  = norm(r_vec);
v2 = dot(v_vec, v_vec);
a  = 1 / (2/R - v2/mu);
if ~isfinite(a) || a <= 0
    rp = NaN;  ra = NaN;  return;
end
h_vec = cross(r_vec, v_vec);
e_vec = cross(v_vec, h_vec)/mu - r_vec/R;
e     = norm(e_vec);
if e >= 1.0
    rp = NaN;  ra = NaN;  return;
end
rp = a*(1-e);
ra = a*(1+e);
end

% -------------------------------------------------------------------------
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
