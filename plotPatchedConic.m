function plotPatchedConic(result, departBody, arrivalBody)
%PLOTPATCHEDCONIC Plot a simple visualization of a patched-conic result
%   plotPatchedConic(result, departBody, arrivalBody)
%
%   This helper plots:
%     - A bar chart of delta-V components (departure, arrival, capture, etc.)
%     - A simple 2D view of the transfer orbit (interplanetary case)

% Plot orbit diagram
figure('Name', sprintf('Patched Conic: %s -> %s', departBody.name, arrivalBody.name), 'NumberTitle', 'off');

if strcmpi(departBody.name, 'Earth') && strcmpi(arrivalBody.name, 'Moon')
	% Earth-Moon patched conic — 3D ecliptic-frame view
	Earth_c    = constants().Earth;
	Moon_c     = constants().Moon;
	i_moon_rad = deg2rad(Moon_c.inclination);   % 5.145° orbital inclination to ecliptic
	i_eq_deg   = 23.45;                          % Earth obliquity to ecliptic
	i_eq_rad   = deg2rad(i_eq_deg);

	% Ascending node placed at +x axis (Ω = 0, simplified model).
	% Moon arrival angle in its orbital plane, departure is 180° opposite.
	theta_arr = deg2rad(result.phaseAngle);
	theta_dep = theta_arr + pi;

	% Moon's full inclined orbit
	th  = linspace(0, 2*pi, 300);
	xMo = Moon_c.a * cos(th);
	yMo = Moon_c.a * sin(th) .* cos(i_moon_rad);
	zMo = Moon_c.a * sin(th) .* sin(i_moon_rad);

	% Moon position at arrival
	xM = Moon_c.a * cos(theta_arr);
	yM = Moon_c.a * sin(theta_arr) * cos(i_moon_rad);
	zM = Moon_c.a * sin(theta_arr) * sin(i_moon_rad);

	% TLI transfer arc (perigee→apogee) in Moon's orbital plane
	r0  = result.details.r0;
	a_t = result.details.transferSemiMajor;
	e_t = (Moon_c.a - r0) / (Moon_c.a + r0);
	nu  = linspace(0, pi, 250);
	r_nu = a_t * (1 - e_t^2) ./ (1 + e_t * cos(nu));
	phi  = theta_dep + nu;
	xTLI = r_nu .* cos(phi);
	yTLI = r_nu .* sin(phi) .* cos(i_moon_rad);
	zTLI = r_nu .* sin(phi) .* sin(i_moon_rad);

	% Parking orbit in Moon's orbital plane
	xPk = r0 * cos(th);
	yPk = r0 * sin(th) .* cos(i_moon_rad);
	zPk = r0 * sin(th) .* sin(i_moon_rad);

	% TLI burn position
	xDep = r0 * cos(theta_dep);
	yDep = r0 * sin(theta_dep) * cos(i_moon_rad);
	zDep = r0 * sin(theta_dep) * sin(i_moon_rad);

	% Reference plane discs
	Rd    = Moon_c.a * 1.12;
	th_d  = linspace(0, 2*pi, 120);
	xEcl  = Rd * cos(th_d);
	yEcl  = Rd * sin(th_d);
	zEcl  = zeros(size(th_d));
	xEq   = Rd * cos(th_d);
	yEq   = Rd * sin(th_d) .* cos(i_eq_rad);
	zEq   = Rd * sin(th_d) .* sin(i_eq_rad);

	% Exaggerated sphere radii for visibility at Moon-distance scale
	R_Ev = Moon_c.a * 0.035;
	R_Mv = Moon_c.a * 0.018;
	[sx, sy, sz] = sphere(28);

	% ---- Draw ----
	patch(xEcl, yEcl, zEcl, [0.92 0.90 0.78], 'FaceAlpha', 0.15, ...
	      'EdgeColor', 'none', 'DisplayName', 'Ecliptic plane');
	hold on;
	patch(xEq, yEq, zEq, [0.55 0.75 1.0], 'FaceAlpha', 0.12, ...
	      'EdgeColor', 'none', 'DisplayName', sprintf('Equatorial plane (ε = %.1f°)', i_eq_deg));

	plot3(xMo, yMo, zMo, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.2, ...
	      'DisplayName', 'Moon orbit');
	plot3(xTLI, yTLI, zTLI, '-', 'Color', [0.6 0.1 0.8], 'LineWidth', 2.5, ...
	      'DisplayName', 'TLI trajectory');
	plot3(xPk, yPk, zPk, '--', 'Color', [0.0 0.75 0.9], 'LineWidth', 1.0, ...
	      'DisplayName', 'Parking orbit');

	surf(sx*R_Ev,        sy*R_Ev,        sz*R_Ev,        'FaceColor', [0.2 0.45 0.75], ...
	     'EdgeColor', 'none', 'DisplayName', 'Earth');
	surf(sx*R_Mv + xM,   sy*R_Mv + yM,   sz*R_Mv + zM,   'FaceColor', [0.72 0.72 0.72], ...
	     'EdgeColor', 'none', 'DisplayName', 'Moon (arrival)');

	plot3(xDep, yDep, zDep, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'c', ...
	      'MarkerEdgeColor', 'k', 'DisplayName', 'TLI burn');
	plot3(xM,   yM,   zM,   'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.1 0.8 0.1], ...
	      'MarkerEdgeColor', 'k', 'DisplayName', 'Lunar arrival');

	% Ecliptic reference frame arrows
	aLen = Moon_c.a * 0.38;
	quiver3(0,0,0, aLen,0,0, 0, 'Color', [0.82 0.72 0.1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
	quiver3(0,0,0, 0,aLen,0, 0, 'Color', [0.82 0.72 0.1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
	quiver3(0,0,0, 0,0,aLen, 0, 'Color', [0.82 0.72 0.1], 'LineWidth', 1.5, 'HandleVisibility', 'off');
	text(aLen*1.04, 0,       0,       '\gamma', 'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.7 0.6 0.05]);
	text(0,       aLen*1.04, 0,       'y_{ecl}',    'FontSize', 8, 'Color', [0.7 0.6 0.05]);
	text(0,       0,       aLen*1.04, 'Ecl. North', 'FontSize', 8, 'Color', [0.7 0.6 0.05]);

	% Ascending node line (Ω = 0, along +x)
	plot3([0 Rd*0.88], [0 0], [0 0], '--', 'Color', [0.55 0.35 0.1], ...
	      'LineWidth', 0.8, 'HandleVisibility', 'off');
	text(Rd*0.90, 0, 0, 'Ascending node', 'FontSize', 7, 'Color', [0.5 0.3 0.1]);

	% Moon orbital inclination arc (in y-z plane: ecliptic y-axis → inclined orbit direction)
	arc_r  = Moon_c.a * 0.22;
	ia     = linspace(0, i_moon_rad, 50);
	plot3(zeros(size(ia)), arc_r*cos(ia), arc_r*sin(ia), '-', ...
	      'Color', [0.5 0.5 0.5], 'LineWidth', 1.1, 'HandleVisibility', 'off');
	text(0, arc_r*cos(i_moon_rad/2)*1.18, arc_r*sin(i_moon_rad/2)*1.18, ...
	     sprintf('i = %.2f°', Moon_c.inclination), 'FontSize', 8, 'Color', [0.4 0.4 0.4]);

	% Earth obliquity arc (smaller radius, same y-z plane)
	ea = linspace(0, i_eq_rad, 50);
	plot3(zeros(size(ea)), arc_r*0.65*cos(ea), arc_r*0.65*sin(ea), '-', ...
	      'Color', [0.3 0.5 0.85], 'LineWidth', 1.1, 'HandleVisibility', 'off');
	text(0, arc_r*0.65*cos(i_eq_rad/2)*1.18, arc_r*0.65*sin(i_eq_rad/2)*1.18, ...
	     sprintf('ε = %.1f°', i_eq_deg), 'FontSize', 8, 'Color', [0.25 0.45 0.75]);

	hold off;
	axis equal;
	grid on;
	xlabel('x (km)  –  Vernal Equinox');
	ylabel('y (km)  –  Ecliptic');
	zlabel('z (km)  –  Ecl. North');
	title('Earth–Moon Transfer (Ecliptic Frame)');
	legend('Location', 'northeastoutside');
	view(32, 22);

elseif isfield(result.details, 'r1') && isfield(result.details, 'r2') && isfield(result.details, 'aTransfer')
	% Interplanetary transfer — work in AU for readable axes
	AU = 149597870.7; % km
	r1 = result.details.r1 / AU;
	r2 = result.details.r2 / AU;
	a  = result.details.aTransfer / AU;
	e  = abs(r2 - r1) / (r2 + r1);

	% Circular planet orbits
	thetaC = linspace(0, 2*pi, 400);
	x1 = r1 * cos(thetaC);
	y1 = r1 * sin(thetaC);
	x2 = r2 * cos(thetaC);
	y2 = r2 * sin(thetaC);

	% Full transfer ellipse (dashed background)
	theta = linspace(0, 2*pi, 400);
	r_ell = a * (1 - e^2) ./ (1 + e*cos(theta));
	x_ell = r_ell .* cos(theta);
	y_ell = r_ell .* sin(theta);

	% Active transfer arc (perihelion -> aphelion, then rotate by phase)
	phase_rad = deg2rad(result.phaseAngle);
	rot = [cos(phase_rad), -sin(phase_rad); sin(phase_rad), cos(phase_rad)];
	transferTheta = linspace(0, pi, 120);
	r_arc = a * (1 - e^2) ./ (1 + e*cos(transferTheta));
	x_arc = r_arc .* cos(transferTheta);
	y_arc = r_arc .* sin(transferTheta);
	rotated = rot * [x_arc; y_arc];

	% Departure and arrival positions (correctly using rotated arc endpoints)
	depPos = rotated(:, 1);
	arrPos = rotated(:, end);

	% --- Draw ---
	plot(0, 0, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [1.0 0.85 0.0], ...
	     'MarkerEdgeColor', [0.8 0.6 0.0], 'DisplayName', 'Sun');
	hold on;
	plot(x1, y1, '-',  'Color', [0.2 0.5 1.0], 'LineWidth', 1.2, 'DisplayName', departBody.name);
	plot(x2, y2, '-',  'Color', [0.9 0.3 0.2], 'LineWidth', 1.2, 'DisplayName', arrivalBody.name);
	plot(x_ell, y_ell, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'DisplayName', 'Transfer ellipse');
	plot(rotated(1,:), rotated(2,:), '-', 'Color', [0.6 0.1 0.8], 'LineWidth', 2.0, 'DisplayName', 'Trajectory');

	% Departure / arrival markers
	plot(depPos(1), depPos(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'c', ...
	     'MarkerEdgeColor', 'k', 'DisplayName', 'Departure');
	plot(arrPos(1), arrPos(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.1 0.8 0.1], ...
	     'MarkerEdgeColor', 'k', 'DisplayName', 'Arrival');

	% Planet name labels near orbits
	text(r1, 0, [' ' departBody.name], 'Color', [0.2 0.5 1.0], 'FontWeight', 'bold', 'FontSize', 9);
	text(r2, 0, [' ' arrivalBody.name], 'Color', [0.9 0.3 0.2], 'FontWeight', 'bold', 'FontSize', 9);

	hold off;
	axis equal;
	title(sprintf('Coplanar Transfer: %s \\rightarrow %s', departBody.name, arrivalBody.name));
	legend('Location', 'best');
	xlabel('AU');
	ylabel('AU');
	grid on;

	% Info annotation: ΔV and TOF
	tofDays = result.tof / 86400;
	infoStr = sprintf('\\DeltaV = %.2f km/s    TOF = %.1f days', result.deltaV, tofDays);
	if isfield(result, 'departureJD')
		depStr  = datestr(result.departureJD - 1721058.5, 'yyyy-mm-dd');
		infoStr = sprintf('%s    Dep: %s', infoStr, depStr);
	end
	text(0.02, 0.02, infoStr, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
	     'FontSize', 8, 'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.7 0.7 0.7]);

	% Body-centric departure and arrival plots
	if isfield(result.details, 'vInfDepart') && isfield(result.details, 'rParkDepart')
		inc_dep = result.details.departureInclination;
		inc_arr = result.details.arrivalInclination;
		figure('Name', sprintf('Body-Centric Views: %s -> %s', departBody.name, arrivalBody.name), 'NumberTitle', 'off', ...
		       'Position', [50 50 1400 680]);

		subplot(1, 2, 1);
		plotBodyCentric3D(departBody, result.details.rParkDepart, result.details.vInfDepart, ...
		    inc_dep, 'hyperbola', sprintf('%s Departure', departBody.name));
		depStr = sprintf('\\DeltaV = %.3f km/s', result.details.dvDeparture);
		if isfield(result, 'departureJD')
			depStr = sprintf('%s\nDep: %s', depStr, datestr(result.departureJD - 1721058.5, 'yyyy-mm-dd'));
		end
		text(0.03, 0.03, depStr, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
		     'FontSize', 8, 'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.7 0.7 0.7]);

		subplot(1, 2, 2);
		plotBodyCentric3D(arrivalBody, result.details.rParkArrive, result.details.vInfArrive, ...
		    inc_arr, 'hyperbola', sprintf('%s Arrival', arrivalBody.name));
		arrStr = sprintf('\\DeltaV = %.3f km/s', result.details.dvArrival);
		if isfield(result, 'departureJD')
			arrJD  = result.departureJD + result.tof / 86400;
			arrStr = sprintf('%s\nArr: %s', arrStr, datestr(arrJD - 1721058.5, 'yyyy-mm-dd'));
		end
		text(0.03, 0.03, arrStr, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
		     'FontSize', 8, 'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.7 0.7 0.7]);
	end
end

% Lunar body-centric plots
if strcmpi(departBody.name, 'Earth') && strcmpi(arrivalBody.name, 'Moon') && ...
        isfield(result.details, 'rParkArrive')
	bodies = constants();
	inc_dep = result.details.departureInclination;
	inc_arr = result.details.arrivalInclination;
	figure('Name', 'Body-Centric Views: Earth -> Moon', 'NumberTitle', 'off', ...
	       'Position', [50 50 1400 680]);

	subplot(1, 2, 1);
	plotBodyCentric3D(bodies.Earth, result.details.r0, result.details.transferSemiMajor, ...
	    inc_dep, 'tli', sprintf('%s Departure (TLI)', bodies.Earth.name));
	view(32+180, 22);  % mirror of arrival view
	depStr = sprintf('\\DeltaV_{TLI} = %.3f km/s', result.details.dvTLI);
	if isfield(result, 'departureJD')
		depStr = sprintf('%s\nDep: %s', depStr, datestr(result.departureJD - 1721058.5, 'yyyy-mm-dd'));
	end
	text(0.03, 0.03, depStr, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
	     'FontSize', 8, 'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.7 0.7 0.7]);

	subplot(1, 2, 2);
	r_apo_moon = 0;
	if isfield(result.details, 'rApoArrive'), r_apo_moon = result.details.rApoArrive; end
	omega_arr = 0;
	if isfield(result.details, 'arrivalArgOfPeriapsis'), omega_arr = result.details.arrivalArgOfPeriapsis; end
	% Build intermediate orbit list for bi-elliptic transfers
	extra_orbits = {};
	if isfield(result.details, 'rBiEllipticApoapsis')
		r_bi  = result.details.rBiEllipticApoapsis;
		r_p   = result.details.rParkArrive;
		alt_bi = r_bi - bodies.Moon.radius;
		extra_orbits{1} = struct('r_peri', r_p, 'r_apo', r_bi, 'inc', 0, ...
		    'label', sprintf('Post-LOI equatorial (apo = %.0f km)', alt_bi));
		extra_orbits{2} = struct('r_peri', r_p, 'r_apo', r_bi, 'inc', inc_arr, ...
		    'label', sprintf('Post-\\Deltai polar (apo = %.0f km)', alt_bi));
	end
	plotBodyCentric3D(bodies.Moon, result.details.rParkArrive, result.details.vInf, ...
	    inc_arr, 'hyperbola', 'Lunar Arrival', r_apo_moon, omega_arr, extra_orbits);
	if isfield(result.details, 'rBiEllipticApoapsis')
		arrStr = sprintf('\\DeltaV_{LOI} = %.3f  \\DeltaV_{\\Deltai} = %.3f  \\DeltaV_{trim} = %.3f km/s', ...
		    result.details.dvLOI, result.details.dvPlaneChange, result.details.dvApoapsisTrim);
	else
		arrStr = sprintf('\\DeltaV_{cap} = %.3f km/s', result.details.dvCapture);
	end
	if isfield(result, 'departureJD')
		arrJD  = result.departureJD + result.tof / 86400;
		arrStr = sprintf('%s\nArr: %s', arrStr, datestr(arrJD - 1721058.5, 'yyyy-mm-dd'));
	end
	text(0.03, 0.03, arrStr, 'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
	     'FontSize', 8, 'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.7 0.7 0.7]);
end

end

% -------------------------------------------------------------------------
function plotBodyCentric3D(body, r_peri, traj_param, inc_deg, traj_type, titleStr, r_apo_orbit, omega_deg, extra_orbits)
%PLOTBODYCENTRIC3D Unified 3D body-centric departure/arrival view.
%   traj_type: 'hyperbola' — draws hyperbolic pass (interplanetary / lunar capture)
%              'tli'       — draws clipped TLI ellipse arc (lunar departure)
%   traj_param: v_inf (km/s) for hyperbola, a_transfer (km) for tli
%   inc_deg: orbit inclination relative to transfer plane (deg)
%   r_apo_orbit: apoapsis radius of capture orbit (km); omit or 0 for circular
%   omega_deg: argument of periapsis (deg); rotates orbit in its own plane
%   extra_orbits: cell array of structs with fields r_peri, r_apo, inc, label
%                 (used for bi-elliptic intermediate orbits)

if nargin < 7 || isempty(r_apo_orbit) || r_apo_orbit <= r_peri
    r_apo_orbit = 0;
end
if nargin < 8 || isempty(omega_deg)
    omega_deg = 0;
end
if nargin < 9 || isempty(extra_orbits)
    extra_orbits = {};
end

% Expand view to contain full orbit when apogee is present
if r_apo_orbit > 0
    lim = max(3 * r_peri, 1.15 * r_apo_orbit);
else
    lim = 3 * r_peri;
end
% Expand further to contain extra orbits (e.g. bi-elliptic intermediate)
for k = 1:numel(extra_orbits)
    lim = max(lim, 1.15 * extra_orbits{k}.r_apo);
end

% --- Trajectory in transfer plane (z = 0) ---
if strcmpi(traj_type, 'hyperbola')
    mu        = body.mu;
    a_h       = mu / traj_param^2;
    e_h       = 1 + r_peri / a_h;
    p_h       = a_h * (e_h^2 - 1);
    theta_inf = acos(-1 / e_h);
    cos_lim   = (p_h / lim - 1) / e_h;
    theta_max = min(acos(max(-1, cos_lim)), 0.99 * theta_inf);
    nu        = linspace(-theta_max, theta_max, 400);
    r_traj    = p_h ./ (1 + e_h * cos(nu));
    trajLabel = 'Hyperbola';
else % tli
    a_transfer = traj_param;
    e          = 1 - r_peri / a_transfer;
    p          = a_transfer * (1 - e^2);
    cos_lim    = (p / lim - 1) / e;
    theta_max  = acos(max(-1, min(1, cos_lim)));
    nu         = linspace(-theta_max, theta_max, 400);
    r_traj     = p ./ (1 + e * cos(nu));
    trajLabel  = 'TLI trajectory';
end
x_traj = r_traj .* cos(nu);
y_traj = r_traj .* sin(nu);
z_traj = zeros(size(nu));

% --- Parking / capture orbit inclined at inc_deg ---
% Ascending node along periapsis direction (+x).
% If apogee is given draw full ellipse, otherwise circular.
inc_rad = deg2rad(inc_deg);
th = linspace(0, 2*pi, 300);
if r_apo_orbit > 0
    a_orb  = (r_peri + r_apo_orbit) / 2;
    e_orb  = (r_apo_orbit - r_peri) / (r_apo_orbit + r_peri);
    r_orb  = a_orb * (1 - e_orb^2) ./ (1 + e_orb * cos(th));
    x_p    = r_orb .* cos(th);
    y_p    = r_orb .* sin(th) * cos(inc_rad);
    z_p    = r_orb .* sin(th) * sin(inc_rad);
    orbitLabel = sprintf('Orbit (i=%.1f°, e=%.3f)', inc_deg, e_orb);
else
    x_p = r_peri * cos(th);
    y_p = r_peri * sin(th) * cos(inc_rad);
    z_p = r_peri * sin(th) * sin(inc_rad);
    orbitLabel = sprintf('Orbit (i=%.1f°)', inc_deg);
end

% --- Argument of periapsis rotation (Rodrigues formula) ---
% Orbit normal for ascending node at +x, inclination i: n = [0; -sin(i); cos(i)]
% Always compute R_om; it is identity when omega_deg == 0.
inc_rad_rot = deg2rad(inc_deg);
n_hat = [0; -sin(inc_rad_rot); cos(inc_rad_rot)];
om    = deg2rad(omega_deg);
nx = n_hat(1); ny = n_hat(2); nz = n_hat(3);
skew_n = [0, -nz,  ny;
          nz,  0, -nx;
         -ny,  nx,  0];
R_om = cos(om)*eye(3) + (1 - cos(om))*(n_hat*n_hat') + sin(om)*skew_n;
if abs(omega_deg) > 0.01
    % Rotate orbit ring
    pts_orb = R_om * [x_p; y_p; z_p];
    x_p = pts_orb(1,:);  y_p = pts_orb(2,:);  z_p = pts_orb(3,:);
    % Rotate trajectory so periapsis aligns with capture orbit periapsis
    pts_traj = R_om * [x_traj; y_traj; z_traj];
    x_traj = pts_traj(1,:);  y_traj = pts_traj(2,:);  z_traj = pts_traj(3,:);
end
peri_rot = R_om * [r_peri; 0; 0];
peri_x = peri_rot(1);  peri_y = peri_rot(2);  peri_z = peri_rot(3);

% --- Reference plane discs ---
th_d = linspace(0, 2*pi, 120);
x_tp = lim * cos(th_d);   y_tp = lim * sin(th_d);   z_tp = zeros(size(th_d));

obl_deg = 0;
if isfield(body, 'obliquity'), obl_deg = body.obliquity; end
obl_rad = deg2rad(obl_deg);
x_eq = lim * cos(th_d);
y_eq = lim * sin(th_d) * cos(obl_rad);
z_eq = lim * sin(th_d) * sin(obl_rad);

% --- Body sphere ---
R_vis    = min(body.radius, lim * 0.33);
[sx,sy,sz] = sphere(28);

% --- Draw ---
patch(x_tp, y_tp, z_tp, [0.92 0.90 0.78], 'FaceAlpha', 0.15, ...
      'EdgeColor', 'none', 'DisplayName', 'Transfer plane');
hold on;
if obl_deg > 0.5
    patch(x_eq, y_eq, z_eq, [0.55 0.75 1.0], 'FaceAlpha', 0.12, ...
          'EdgeColor', 'none', 'DisplayName', sprintf('Equatorial plane (ε=%.1f°)', obl_deg));
end
surf(sx*R_vis, sy*R_vis, sz*R_vis, 'FaceColor', bodyColor(body.name), ...
     'EdgeColor', 'none', 'DisplayName', body.name);
plot3(x_p, y_p, z_p, '--', 'Color', [0.0 0.75 0.9], 'LineWidth', 1.2, ...
      'DisplayName', orbitLabel);
plot3(x_traj, y_traj, z_traj, '-', 'Color', [0.6 0.1 0.8], 'LineWidth', 2, ...
      'DisplayName', trajLabel);
plot3(peri_x, peri_y, peri_z, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'c', ...
      'MarkerEdgeColor', 'k', 'DisplayName', 'Periapsis (Burn 1 & 3)');

% --- Bi-elliptic intermediate orbits ---
% Each extra orbit is rotated by the same R_om so all periapses coincide.
% Orbit 1 (equatorial, high apo): drawn in the transfer plane (inc=0) then rotated.
% Orbit 2 (polar,      high apo): drawn with full inc then rotated — same shape,
%   different plane → visually distinct, showing the plane change.
ex_colors = {[0.95 0.55 0.05], [0.15 0.82 0.25]};   % orange, lime-green
ex_styles  = {'--', ':'};
for k = 1:numel(extra_orbits)
    eo    = extra_orbits{k};
    th_e  = linspace(0, 2*pi, 300);
    a_e   = (eo.r_peri + eo.r_apo) / 2;
    e_e   = (eo.r_apo - eo.r_peri) / (eo.r_apo + eo.r_peri);
    r_e   = a_e * (1 - e_e^2) ./ (1 + e_e * cos(th_e));
    inc_e = deg2rad(eo.inc);
    xe = r_e .* cos(th_e);
    ye = r_e .* sin(th_e) * cos(inc_e);
    ze = r_e .* sin(th_e) * sin(inc_e);
    pts_e = R_om * [xe; ye; ze];
    xe = pts_e(1,:);  ye = pts_e(2,:);  ze = pts_e(3,:);
    cidx = min(k, numel(ex_colors));
    sidx = min(k, numel(ex_styles));
    plot3(xe, ye, ze, ex_styles{sidx}, 'Color', ex_colors{cidx}, ...
          'LineWidth', 1.2, 'DisplayName', eo.label);
end
% Burn 2 marker: plane change occurs at apolune of the intermediate orbit.
% After R_om the apolune of the equatorial intermediate orbit (nu = π → −r_bi_apo, 0, 0)
% maps to (0, 0, −r_bi_apo) regardless of inclination — directly above the south pole.
if numel(extra_orbits) >= 1
    apo_pt = R_om * [-extra_orbits{1}.r_apo; 0; 0];
    plot3(apo_pt(1), apo_pt(2), apo_pt(3), 's', 'MarkerSize', 9, ...
          'MarkerFaceColor', [1.0 0.85 0.0], 'MarkerEdgeColor', 'k', ...
          'DisplayName', sprintf('Burn 2 – plane change (%.0f km alt)', ...
          extra_orbits{1}.r_apo - body.radius));
end

% North/South pole labels when orbit is significantly polar and ω is set
if inc_deg > 45 && abs(omega_deg) > 0.01
    text(0, 0,  lim*0.55, 'North Pole', 'FontSize', 8, 'Color', [0.2 0.6 0.2], ...
         'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
    text(0, 0, -lim*0.55, 'South Pole', 'FontSize', 8, 'Color', [0.55 0.35 0.1], ...
         'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
    % Landing zone marker at south pole surface
    text(0, 0, -body.radius*1.15, '\downarrow Landing Zone', 'FontSize', 8, ...
         'Color', [0.85 0.25 0.1], 'HorizontalAlignment', 'center', ...
         'FontWeight', 'bold', 'HandleVisibility', 'off');
end

% Reference frame arrows
aLen = lim * 0.42;
quiver3(0,0,0, aLen,0,0, 0, 'Color',[0.82 0.72 0.1],'LineWidth',1.5,'HandleVisibility','off');
quiver3(0,0,0, 0,aLen,0, 0, 'Color',[0.82 0.72 0.1],'LineWidth',1.5,'HandleVisibility','off');
quiver3(0,0,0, 0,0,aLen, 0, 'Color',[0.82 0.72 0.1],'LineWidth',1.5,'HandleVisibility','off');
text(aLen*1.05, 0,       0,       'Periapsis', 'FontSize',7,'Color',[0.7 0.6 0.05]);
text(0,       aLen*1.05, 0,       'y',         'FontSize',8,'Color',[0.7 0.6 0.05]);
text(0,       0,       aLen*1.05, 'Normal',    'FontSize',8,'Color',[0.7 0.6 0.05]);

% Orbit inclination arc (y-z plane: transfer-plane y-axis → inclined orbit)
arc_r = lim * 0.22;
if inc_deg > 0.5
    ia = linspace(0, inc_rad, 50);
    plot3(zeros(size(ia)), arc_r*cos(ia), arc_r*sin(ia), '-', ...
          'Color',[0.5 0.5 0.5],'LineWidth',1.1,'HandleVisibility','off');
    text(0, arc_r*cos(inc_rad/2)*1.18, arc_r*sin(inc_rad/2)*1.18, ...
         sprintf('i=%.1f°', inc_deg), 'FontSize',8,'Color',[0.4 0.4 0.4]);
end

% Obliquity arc
if obl_deg > 0.5
    ea = linspace(0, obl_rad, 50);
    plot3(zeros(size(ea)), arc_r*0.65*cos(ea), arc_r*0.65*sin(ea), '-', ...
          'Color',[0.3 0.5 0.85],'LineWidth',1.1,'HandleVisibility','off');
    text(0, arc_r*0.65*cos(obl_rad/2)*1.18, arc_r*0.65*sin(obl_rad/2)*1.18, ...
         sprintf('ε=%.1f°', obl_deg), 'FontSize',8,'Color',[0.25 0.45 0.75]);
end

hold off;
axis equal;
grid on;
xlim([-lim lim]); ylim([-lim lim]); zlim([-lim lim]);
xlabel('km'); ylabel('km'); zlabel('km');
legend('Location', 'northeast', 'FontSize', 7);
title(titleStr);
view(32, 22);
end

% -------------------------------------------------------------------------
function col = bodyColor(name)
switch lower(name)
    case 'earth',   col = [0.20 0.45 0.75];
    case 'mars',    col = [0.72 0.28 0.18];
    case 'moon',    col = [0.72 0.72 0.72];
    case 'venus',   col = [0.85 0.75 0.35];
    case 'jupiter', col = [0.75 0.60 0.45];
    case 'ceres',   col = [0.55 0.52 0.50];
    case 'vesta',   col = [0.60 0.58 0.52];
    otherwise,      col = [0.50 0.55 0.60];
end
end
