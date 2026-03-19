function plotFlybySequence(result, bodies)
%PLOTFLYBYSEQUENCE  3D heliocentric visualization of a gravity-assist sequence.
%
%   plotFlybySequence(result, bodies)
%
%   Inputs:
%     result - output struct from flybySequence()
%     bodies - 1×N cell array of body structs (same order passed to flybySequence)
%
%   Produces a single 3D figure showing:
%     - Ecliptic reference plane (transparent gray disk)
%     - Sun at the origin
%     - Full orbit ellipse for each body in the sequence
%     - Active Lambert arc for each leg (colored by leg index)
%     - Body positions at each encounter (filled spheres)
%     - Departure body position at arrival (dashed marker)
%     - Annotation box with ΔV budget and flyby summary

N     = numel(bodies);
nLegs = N - 1;
nFB   = N - 2;
muSun = constants().Sun.mu;

legColors = [0.20 0.50 1.00;   % blue
             0.15 0.75 0.30;   % green
             0.95 0.55 0.10;   % orange
             0.80 0.15 0.15;   % red
             0.60 0.20 0.80;   % purple
             0.10 0.80 0.80];  % cyan
if nLegs > size(legColors,1)
    legColors = repmat(legColors, ceil(nLegs/size(legColors,1)), 1);
end

fig = figure('Name', sprintf('Gravity Assist: %s', sequenceName(bodies)), ...
    'NumberTitle', 'off', 'Color', [0.06 0.06 0.10]);
ax = axes('Parent', fig, 'Color', [0.06 0.06 0.10], ...
    'XColor', [0.7 0.7 0.7], 'YColor', [0.7 0.7 0.7], 'ZColor', [0.7 0.7 0.7]);
hold(ax, 'on');  axis(ax, 'equal');  grid(ax, 'on');
ax.GridColor = [0.3 0.3 0.3];
xlabel(ax, 'X_{ecl} (km)');
ylabel(ax, 'Y_{ecl} (km)');
zlabel(ax, 'Z_{ecl} (km)');
title(ax, sprintf('Gravity-Assist Trajectory: %s', sequenceName(bodies)), ...
    'Color', [0.9 0.9 0.9]);

% --- Ecliptic reference plane ---
rMax = 0;
for i = 1:N
    if isfield(bodies{i}, 'a'), rMax = max(rMax, bodies{i}.a); end
end
[xe, ye] = meshgrid(linspace(-rMax, rMax, 3));
surf(ax, xe, ye, zeros(size(xe)), 'FaceColor', [0.4 0.4 0.5], ...
    'FaceAlpha', 0.08, 'EdgeColor', 'none');

% --- Sun ---
R_sun = rMax * 0.015;
[xs, ys, zs] = sphere(20);
surf(ax, xs*R_sun, ys*R_sun, zs*R_sun, 'FaceColor', [1.0 0.9 0.2], ...
    'EdgeColor', 'none', 'FaceAlpha', 1.0);

% --- Planet orbit ellipses ---
for i = 1:N
    [xo, yo, zo] = planetOrbit3D(bodies{i});
    col = bodyColor(bodies{i}.name);
    plot3(ax, xo, yo, zo, '-', 'Color', [col 0.45], 'LineWidth', 1.0);
end

% --- Lambert arcs and encounter markers ---
legHandles = gobjects(nLegs, 1);
for i = 1:nLegs
    leg = result.legs(i);
    col = legColors(i,:);

    % Active Lambert arc
    [xa, ya, za] = transferArc3D(leg.r1_vec, leg.v1_transfer, leg.r2_vec, muSun);
    legHandles(i) = plot3(ax, xa, ya, za, '-', 'Color', col, 'LineWidth', 2.0);

    % Full transfer ellipse (dashed)
    [xf, yf, zf] = fullOrbitTrace(leg.r1_vec, leg.v1_transfer, muSun);
    plot3(ax, xf, yf, zf, '--', 'Color', [col 0.25], 'LineWidth', 0.8);
end

% --- Body positions at encounter ---
% Departure body at departure
r1 = result.legs(1).r1_vec;
plot3(ax, r1(1), r1(2), r1(3), 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', bodyColor(bodies{1}.name), ...
    'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
text(ax, r1(1), r1(2), r1(3), sprintf('  %s dep', bodies{1}.name), ...
    'Color', 'w', 'FontSize', 8);

% Flyby bodies
for j = 1:nFB
    r_fb = result.legs(j).r2_vec;
    col  = bodyColor(bodies{j+1}.name);
    plot3(ax, r_fb(1), r_fb(2), r_fb(3), 'd', 'MarkerSize', 11, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
    fb = result.flybys(j);
    lbl = sprintf('  %s flyby\n  v∞=%.2f km/s  δ=%.1f°', ...
        bodies{j+1}.name, fb.v_inf_in, fb.deflection);
    text(ax, r_fb(1), r_fb(2), r_fb(3), lbl, 'Color', [0.9 0.9 0.7], 'FontSize', 7);
end

% Arrival body at arrival
rN = result.legs(nLegs).r2_vec;
plot3(ax, rN(1), rN(2), rN(3), 's', 'MarkerSize', 11, ...
    'MarkerFaceColor', bodyColor(bodies{N}.name), ...
    'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
text(ax, rN(1), rN(2), rN(3), sprintf('  %s arr', bodies{N}.name), ...
    'Color', 'w', 'FontSize', 8);

% Departure body position at arrival time (for context)
jd_arr = result.legs(nLegs).arriveJD;
try
    [r_dep_arr, ~] = orbitalState(bodies{1}, jd_arr);
    col = bodyColor(bodies{1}.name);
    plot3(ax, r_dep_arr(1), r_dep_arr(2), r_dep_arr(3), 'x', ...
        'MarkerSize', 8, 'Color', [col 0.5], 'LineWidth', 1.5);
catch
end

% --- Legend ---
legLabels = cell(nLegs, 1);
for i = 1:nLegs
    legLabels{i} = sprintf('Leg %d: %s→%s  (%.0f d)', i, ...
        bodies{i}.name, bodies{i+1}.name, result.legs(i).tofDays);
end
legend(ax, legHandles, legLabels, 'Location', 'best', ...
    'TextColor', [0.85 0.85 0.85], 'Color', [0.12 0.12 0.16], ...
    'EdgeColor', [0.4 0.4 0.4], 'FontSize', 8);

% --- Info annotation ---
d = result.details;
jd0     = d.departureJD;
jdArriv = jd0 + result.tof;
dep_str = jd2datestr(jd0);
arr_str = jd2datestr(jdArriv);

lines{1} = sprintf('Departure:  %s  (%s)', bodies{1}.name, dep_str);
lines{2} = sprintf('Arrival:    %s  (%s)', bodies{N}.name, arr_str);
lines{3} = sprintf('Total TOF:  %.1f days', result.tof);
lines{4} = '';
lines{5} = sprintf('C3:             %.2f km²/s²', d.C3);
lines{6} = sprintf('v∞ dep:         %.3f km/s',  d.vInfDepart);
lines{7} = sprintf('v∞ arr:         %.3f km/s',  d.vInfArrive);
lines{8} = '';
lines{9} = sprintf('ΔV departure:   %.3f km/s',  d.dvDeparture);
if nFB > 0
    lines{10} = sprintf('ΔV powered GA:  %.3f km/s', d.dvPoweredFlybys);
    lines{11} = sprintf('ΔV arrival:     %.3f km/s', d.dvArrival);
    lines{12} = sprintf('ΔV TCM:         %.0f m/s',  d.dvTCM*1000);
    lines{13} = sprintf('ΔV total:       %.3f km/s', result.deltaV);
else
    lines{10} = sprintf('ΔV arrival:     %.3f km/s', d.dvArrival);
    lines{11} = sprintf('ΔV TCM:         %.0f m/s',  d.dvTCM*1000);
    lines{12} = sprintf('ΔV total:       %.3f km/s', result.deltaV);
    lines{13} = '';
end
annotation(fig, 'textbox', [0.01 0.01 0.30 0.40], ...
    'String', lines(~cellfun(@(s) isempty(s) && false, lines)), ...
    'FontSize', 8, 'FontName', 'Courier', ...
    'BackgroundColor', [0.10 0.10 0.15], 'EdgeColor', [0.40 0.40 0.55], ...
    'Color', [0.85 0.85 0.90], 'FitBoxToText', true, 'Interpreter', 'none');

view(ax, 25, 28);
end

% =========================================================================
%  Local helpers
% =========================================================================

function name = sequenceName(bodies)
names = cellfun(@(b) b.name, bodies, 'UniformOutput', false);
name  = strjoin(names, ' → ');
end

function s = jd2datestr(jd)
% Convert Julian Date to 'YYYY-Mon-DD' string (no toolbox needed)
jd0 = jd - 1721058.5;   % offset to MATLAB serial date base (~Jan 0, 0000)
try
    s = datestr(jd0, 'yyyy-mmm-dd');
catch
    s = sprintf('JD %.1f', jd);
end
end

% -------------------------------------------------------------------------
function [x, y, z] = planetOrbit3D(body)
if ~isfield(body,'Omega') || ~isfield(body,'omega_peri') || ~isfield(body,'M0')
    th = linspace(0, 2*pi, 360);
    x  = body.a * cos(th);  y = body.a * sin(th);  z = zeros(size(th));
    return;
end
a   = body.a;  e = body.e;
Om  = deg2rad(body.Omega);
om  = deg2rad(body.omega_peri);
inc = deg2rad(body.inclination);
cOm = cos(Om);  sOm = sin(Om);
com = cos(om);  som = sin(om);
ci  = cos(inc); si  = sin(inc);
R = [ cOm*com - sOm*som*ci,  -cOm*som - sOm*com*ci,  sOm*si;
      sOm*com + cOm*som*ci,  -sOm*som + cOm*com*ci, -cOm*si;
      si*som,                  si*com,                 ci    ];
nu    = linspace(0, 2*pi, 360);
p     = a * (1 - e^2);
r     = p ./ (1 + e*cos(nu));
r_pqw = [r.*cos(nu); r.*sin(nu); zeros(1,360)];
r_ecl = R * r_pqw;
x = r_ecl(1,:);  y = r_ecl(2,:);  z = r_ecl(3,:);
end

% -------------------------------------------------------------------------
function [x, y, z] = fullOrbitTrace(r_vec, v_vec, mu)
r_vec = r_vec(:);  v_vec = v_vec(:);
R1    = norm(r_vec);
h_vec = cross(r_vec, v_vec);
e_vec = cross(v_vec, h_vec)/mu - r_vec/R1;
e     = norm(e_vec);
if e > 1e-8, P_hat = e_vec/e; else, P_hat = r_vec/R1; end
W_hat = h_vec / norm(h_vec);
Q_hat = cross(W_hat, P_hat);
a_t   = 1 / (2/R1 - dot(v_vec,v_vec)/mu);
p_t   = a_t * (1 - e^2);
nu    = linspace(0, 2*pi, 360);
r     = p_t ./ (1 + e*cos(nu));
x = r .* (P_hat(1)*cos(nu) + Q_hat(1)*sin(nu));
y = r .* (P_hat(2)*cos(nu) + Q_hat(2)*sin(nu));
z = r .* (P_hat(3)*cos(nu) + Q_hat(3)*sin(nu));
end

% -------------------------------------------------------------------------
function [x, y, z] = transferArc3D(r1_vec, v1_vec, r2_vec, mu)
r1_vec = r1_vec(:);  v1_vec = v1_vec(:);  r2_vec = r2_vec(:);
R1    = norm(r1_vec);
h_vec = cross(r1_vec, v1_vec);
e_vec = cross(v1_vec, h_vec)/mu - r1_vec/R1;
e     = norm(e_vec);
if e > 1e-8, P_hat = e_vec/e; else, P_hat = r1_vec/R1; end
W_hat = h_vec / norm(h_vec);
Q_hat = cross(W_hat, P_hat);
a_t   = 1 / (2/R1 - dot(v1_vec,v1_vec)/mu);
p_t   = a_t * (1 - e^2);
nu1   = atan2(dot(r1_vec, Q_hat), dot(r1_vec, P_hat));
nu2   = atan2(dot(r2_vec, Q_hat), dot(r2_vec, P_hat));
if nu2 <= nu1, nu2 = nu2 + 2*pi; end
nu    = linspace(nu1, nu2, 300);
r     = p_t ./ (1 + e*cos(nu));
x = r .* (P_hat(1)*cos(nu) + Q_hat(1)*sin(nu));
y = r .* (P_hat(2)*cos(nu) + Q_hat(2)*sin(nu));
z = r .* (P_hat(3)*cos(nu) + Q_hat(3)*sin(nu));
end

% -------------------------------------------------------------------------
function col = bodyColor(name)
switch lower(name)
    case 'mercury',  col = [0.60 0.58 0.55];
    case 'earth',    col = [0.20 0.45 0.75];
    case 'mars',     col = [0.72 0.28 0.18];
    case 'moon',     col = [0.72 0.72 0.72];
    case 'venus',    col = [0.85 0.75 0.35];
    case 'jupiter',  col = [0.75 0.60 0.45];
    case 'saturn',   col = [0.85 0.80 0.55];
    case 'uranus',   col = [0.55 0.85 0.85];
    case 'neptune',  col = [0.25 0.40 0.85];
    case 'pluto',    col = [0.70 0.60 0.55];
    otherwise,       col = [0.70 0.70 0.70];
end
end
