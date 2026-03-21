% example_100au.m
% Mission design: reach 100 AU in 10 years.
%
% Architecture:
%   1. Earth departure (C3 launch)
%   2. Jupiter gravity assist — max retrograde turn, spacecraft dives toward Sun
%   3. Solar Oberth maneuver at perihelion (floor = Parker Solar Probe limit)
%   4. Hyperbolic escape to 100 AU
%
% Thermal constraint: perihelion >= 9.86 solar radii (PSP planned minimum).
% Direction of escape is unconstrained — goal is purely to minimize transit time.

bodies = constants();
AU    = bodies.Constants.AU;
muSun = bodies.Sun.mu;

R_sun      = bodies.Sun.radius;          % km
r_peri_min = 9.86 * R_sun;              % km, PSP thermal floor

fprintf('\n======================================================\n');
fprintf('  100 AU in 10 Years: Jupiter + Solar Oberth Design\n');
fprintf('======================================================\n\n');
fprintf('Thermal floor: %.2f R_sun = %.4f AU = %.0f km\n\n', ...
    r_peri_min/R_sun, r_peri_min/AU, r_peri_min);

%% ---- Parameters --------------------------------------------------------
dv_oberth = 6.0;    % km/s  Oberth kick stage (chemical, Isp = 450 s)
isp_chem  = 450;    % s

r_target    = 100 * AU;           % km
t_goal_yr   = 10;                 % years

%% ---- Design scan: launch date vs Jupiter flyby date --------------------
fprintf('Scanning launch/flyby grid...\n');

jd_launch_start = julianDate(2025, 1, 1);
jd_launch_end   = julianDate(2033, 6, 1);
jd_jup_start    = julianDate(2026, 1, 1);
jd_jup_end      = julianDate(2037, 6, 1);

nL = 60;  nJ = 60;
jd_launches = linspace(jd_launch_start, jd_launch_end, nL);
jd_jups     = linspace(jd_jup_start,    jd_jup_end,    nJ);

C3_grid     = nan(nL, nJ);
rperi_grid  = nan(nL, nJ);   % AU
vinf_grid   = nan(nL, nJ);   % km/s post-Oberth
ttot_grid   = nan(nL, nJ);   % years, launch to 100 AU

r_jup_body = bodies.Jupiter.radius;
mu_jup     = bodies.Jupiter.mu;

for iL = 1:nL
    jdL = jd_launches(iL);
    [r_earth, v_earth] = orbitalState(bodies.Earth, jdL);

    for iJ = 1:nJ
        jdJ    = jd_jups(iJ);
        tof_d  = jdJ - jdL;
        if tof_d < 200 || tof_d > 3*365.25, continue; end

        [r_jup, v_jup] = orbitalState(bodies.Jupiter, jdJ);
        try
            [v1, v2] = lambertSolver(r_earth, r_jup, tof_d*86400, muSun);
        catch, continue; end

        C3 = norm(v1 - v_earth)^2;
        if C3 > 200, continue; end

        % Jupiter arrival v_inf
        vinf_arr_vec = v2 - v_jup;
        vinf_arr     = norm(vinf_arr_vec);

        % Max retrograde turn at Jupiter (periapsis = 1.05 R_Jup)
        rp_jup  = 1.05 * r_jup_body;
        delta   = 2 * asin(1 / (1 + rp_jup * vinf_arr^2 / mu_jup));
        vhat_in = vinf_arr_vec / vinf_arr;
        tgt     = -v_jup / norm(v_jup);           % retrograde = toward Sun
        ang     = acos(max(-1, min(1, dot(vhat_in, tgt))));
        bend    = min(delta, ang);
        perp    = tgt - vhat_in * dot(tgt, vhat_in);
        if norm(perp) > 1e-10
            vhat_out = vhat_in*cos(bend) + (perp/norm(perp))*sin(bend);
        else
            vhat_out = vhat_in;
        end

        v_post = v_jup + vhat_out * vinf_arr;
        r_jd   = norm(r_jup);
        eps    = 0.5*norm(v_post)^2 - muSun/r_jd;

        % Post-flyby perihelion
        if eps >= 0
            r_peri = r_peri_min;
        else
            a_pf  = -muSun / (2*eps);
            h_pf  = norm(cross(r_jup, v_post));
            e_pf  = sqrt(max(0, 1 - h_pf^2/(muSun*a_pf)));
            r_peri = max(a_pf*(1-e_pf), r_peri_min);
        end

        % Speed at perihelion before/after Oberth
        if eps < 0
            vpre = sqrt(muSun*(2/r_peri - 1/(-muSun/(2*eps))));
        else
            vpre = sqrt(norm(v_post)^2 - 2*muSun/r_jd + 2*muSun/r_peri);
        end
        vesc  = sqrt(2*muSun/r_peri);
        if (vpre + dv_oberth) < vesc, continue; end
        vinf_post = sqrt((vpre + dv_oberth)^2 - vesc^2);

        % TOF: leg 1 (Lambert) + leg 2 (Jupiter→perihelion) + leg 3 (perihelion→100 AU)
        tof1 = tof_d * 86400;
        tof2 = 0;
        if eps < 0
            tof2 = pi * sqrt((-muSun/(2*eps))^3 / muSun);   % ≈ half period
        end
        r_arr = linspace(r_peri, r_target, 5000);
        v_arr = sqrt(max(0, 2*(0.5*vinf_post^2 + muSun./r_arr)));
        tof3  = sum(gradient(r_arr) ./ v_arr);

        C3_grid(iL,iJ)    = C3;
        rperi_grid(iL,iJ) = r_peri / AU;
        vinf_grid(iL,iJ)  = vinf_post;
        ttot_grid(iL,iJ)  = (tof1+tof2+tof3) / (365.25*86400);
    end
end
fprintf('Done.\n\n');

%% ---- Best trajectory ---------------------------------------------------
score = ttot_grid;
score(C3_grid > 150 | isnan(score)) = inf;
[~, idx] = min(score(:));
[iL_best, iJ_best] = ind2sub([nL, nJ], idx);

jd_launch_best = jd_launches(iL_best);
jd_jup_best    = jd_jups(iJ_best);
tof_leg1_days  = jd_jup_best - jd_launch_best;

fprintf('Best trajectory (dv_oberth = %.1f km/s):\n', dv_oberth);
fprintf('  Launch:           %s\n', datestr(jd_launch_best - 1721058.5,'dd mmm yyyy'));
fprintf('  Jupiter flyby:    %s  (TOF = %.0f days)\n', ...
    datestr(jd_jup_best - 1721058.5,'dd mmm yyyy'), tof_leg1_days);
fprintf('  C3:               %.1f km^2/s^2\n',   C3_grid(iL_best,iJ_best));
fprintf('  Perihelion:       %.4f AU (%.1f R_sun)\n', ...
    rperi_grid(iL_best,iJ_best), rperi_grid(iL_best,iJ_best)*AU/R_sun);
fprintf('  v_inf post-Oberth: %.2f km/s\n',      vinf_grid(iL_best,iJ_best));
fprintf('  Time to 100 AU:   %.2f years\n\n',    ttot_grid(iL_best,iJ_best));

%% ---- Detailed computation of best trajectory ---------------------------
[r_earth_dep, v_earth_dep] = orbitalState(bodies.Earth,   jd_launch_best);
[r_jup_arr,   v_jup_arr]   = orbitalState(bodies.Jupiter, jd_jup_best);
[v1_best, v2_best] = lambertSolver(r_earth_dep, r_jup_arr, tof_leg1_days*86400, muSun);

vinf_arr_vec = v2_best - v_jup_arr;
vinf_arr_mag = norm(vinf_arr_vec);

% Jupiter flyby
rp_fly     = 1.05 * r_jup_body;
delta_fly  = 2 * asin(1 / (1 + rp_fly * vinf_arr_mag^2 / mu_jup));
vhat_in    = vinf_arr_vec / vinf_arr_mag;
vjup_hat   = v_jup_arr / norm(v_jup_arr);
tgt        = -vjup_hat;
ang        = acos(max(-1, min(1, dot(vhat_in, tgt))));
bend       = min(delta_fly, ang);
perp       = tgt - vhat_in * dot(tgt, vhat_in);
if norm(perp) > 1e-10
    vhat_out_best = vhat_in*cos(bend) + (perp/norm(perp))*sin(bend);
else
    vhat_out_best = vhat_in;
end
v_helio_post_jup = v_jup_arr + vhat_out_best * vinf_arr_mag;

% Post-flyby orbit elements
r_jup_dist = norm(r_jup_arr);
eps_post   = 0.5*norm(v_helio_post_jup)^2 - muSun/r_jup_dist;
a_post     = -muSun/(2*eps_post);
h_post_vec = cross(r_jup_arr, v_helio_post_jup);
p_post     = norm(h_post_vec)^2 / muSun;
e_post     = sqrt(max(0, 1 - p_post/a_post));
r_peri_post   = a_post*(1-e_post);
r_peri_actual = max(r_peri_post, r_peri_min);

v_peri     = sqrt(muSun*(2/r_peri_actual - 1/a_post));
v_esc_peri = sqrt(2*muSun/r_peri_actual);

fprintf('Jupiter flyby details:\n');
fprintf('  v_inf at Jupiter:  %.3f km/s\n', vinf_arr_mag);
fprintf('  Turn angle:        %.1f deg\n', delta_fly*180/pi);
fprintf('  Periapsis altitude: %.0f km (%.2f R_Jup)\n', rp_fly-r_jup_body, rp_fly/r_jup_body);

fprintf('\nPost-flyby orbit:\n');
fprintf('  a = %.4f AU,  e = %.5f\n', a_post/AU, e_post);
fprintf('  Perihelion: %.4f AU (%.2f R_sun)\n', r_peri_actual/AU, r_peri_actual/R_sun);
fprintf('  v_peri pre-burn:  %.2f km/s\n', v_peri);
fprintf('  v_esc at peri:    %.2f km/s\n', v_esc_peri);

%% ---- Oberth ΔV trade ---------------------------------------------------
tof1_best = tof_leg1_days * 86400;
tof2_best = pi * sqrt(a_post^3 / muSun);   % Jupiter→perihelion (≈ half period)

dv_sweep = linspace(0, 15, 500);
vinf_sweep = real(sqrt(max(0, (v_peri + dv_sweep).^2 - v_esc_peri^2)));
ttot_sweep = nan(size(dv_sweep));
for i = 1:numel(dv_sweep)
    if vinf_sweep(i) < 0.1, continue; end
    r_a = linspace(r_peri_actual, r_target, 5000);
    v_a = sqrt(max(0, 2*(0.5*vinf_sweep(i)^2 + muSun./r_a)));
    ttot_sweep(i) = (tof1_best + tof2_best + sum(gradient(r_a)./v_a)) / (365.25*86400);
end

% ΔV for 10-year goal
dv_for_10yr = NaN;
if any(ttot_sweep <= t_goal_yr & ~isnan(ttot_sweep))
    dv_for_10yr = dv_sweep(find(ttot_sweep <= t_goal_yr, 1));
end

fprintf('\nOberth trade at best trajectory perihelion (%.4f AU):\n', r_peri_actual/AU);
fprintf('  %-10s  %-20s  %-16s  %-14s\n', 'dv (km/s)', 'v_inf post (km/s)', 'Time to 100 AU', 'Prop frac');
fprintf('  %s\n', repmat('-', 1, 63));
for dv_s = [1 2 3 5 7 10 12 15]
    vi = sqrt(max(0, (v_peri+dv_s)^2 - v_esc_peri^2));
    if vi < 0.1
        fprintf('  %-10.0f  (below escape)\n', dv_s); continue
    end
    r_a = linspace(r_peri_actual, r_target, 5000);
    v_a = sqrt(max(0, 2*(0.5*vi^2 + muSun./r_a)));
    tt  = (tof1_best + tof2_best + sum(gradient(r_a)./v_a)) / (365.25*86400);
    pf  = 1 - exp(-dv_s / (isp_chem * 9.80665e-3));
    fprintf('  %-10.0f  %-20.2f  %-16.1f  %.1f%%\n', dv_s, vi, tt, 100*pf);
end
if ~isnan(dv_for_10yr)
    fprintf('\n  => %.1f km/s Oberth burn achieves 10-year goal (100 AU)\n\n', dv_for_10yr);
else
    fprintf('\n  => 10-year goal requires >15 km/s Oberth burn at this perihelion\n\n');
end

%% ---- Comparison: no Oberth, direct Jupiter escape ----------------------
fprintf('Context comparison:\n');
fprintf('  New Horizons (direct, v_inf~14 km/s):  ~36 years to 100 AU\n');
fprintf('  Voyager 1    (direct, v_inf~17 km/s):  ~35 years to 100 AU\n');
for vi_ref = [14, 17, vinf_grid(iL_best,iJ_best)]
    r_a = linspace(r_peri_actual, r_target, 5000);
    v_a = sqrt(max(0, 2*(0.5*vi_ref^2 + muSun./r_a)));
    tt  = sum(gradient(r_a)./v_a) / (365.25*86400);
    fprintf('  This mission, v_inf = %.1f km/s (coast from peri): %.1f years to 100 AU\n', vi_ref, tt);
end
fprintf('\n');

%% ---- Pre-compute arc data (shared across plots) -------------------------
e_vec_pl  = cross(v_helio_post_jup, h_post_vec)/muSun - r_jup_arr/norm(r_jup_arr);
e_hat_pl  = e_vec_pl / norm(e_vec_pl);
h_hat_pl  = h_post_vec / norm(h_post_vec);
q_hat_pl  = cross(h_hat_pl, e_hat_pl);
R2d       = [e_hat_pl(1), q_hat_pl(1); e_hat_pl(2), q_hat_pl(2)];

r_jup_hat = r_jup_arr / norm(r_jup_arr);
nu_jup    = atan2(dot(q_hat_pl, r_jup_hat), dot(e_hat_pl, r_jup_hat));
if nu_jup >= 0, nu_end_in = 2*pi; else, nu_end_in = 0; end
nu_in  = linspace(nu_jup, nu_end_in, 1200);
p_in   = a_post * (1 - e_post^2);
r_in   = p_in ./ (1 + e_post * cos(nu_in));
xy_in  = R2d * [r_in .* cos(nu_in); r_in .* sin(nu_in)] / AU;

r_peri_vec  = e_hat_pl * r_peri_actual;
v_peri_dir  = cross(h_hat_pl, e_hat_pl);
v_peri_dir  = v_peri_dir / norm(v_peri_dir);
vinf_actual = sqrt(max(0, (v_peri + dv_oberth)^2 - v_esc_peri^2));
e_esc       = 1 + r_peri_actual * vinf_actual^2 / muSun;
p_esc       = r_peri_actual * (1 + e_esc);
nu_max_esc  = acos(-1/e_esc) * 0.88;
nu_esc      = linspace(0, nu_max_esc, 1200);
r_esc       = p_esc ./ (1 + e_esc * cos(nu_esc));
xy_esc      = R2d * [r_esc .* cos(nu_esc); r_esc .* sin(nu_esc)] / AU;

%% ---- Plots -------------------------------------------------------------
bgCol = [0.06 0.06 0.10];
txtCol= [0.85 0.85 0.85];
axCol = [0.68 0.68 0.68];
gridC = [0.20 0.20 0.26];

th       = linspace(0, 2*pi, 360);
r_sun_AU = bodies.Sun.radius / AU;

%% Plot 1: Oberth ΔV trade ------------------------------------------------
fig1 = figure('Name','Oberth Trade','NumberTitle','off', ...
    'Color',bgCol,'Position',[80 80 720 440]);
ax1 = axes('Parent',fig1,'Color',bgCol,'XColor',axCol,'YColor',axCol, ...
    'GridColor',gridC,'Box','on');
hold(ax1,'on');  grid(ax1,'on');

yyaxis(ax1,'left');
plot(ax1, dv_sweep, ttot_sweep, '-', 'Color',[0.35 0.75 1.00], 'LineWidth',2.5);
yline(ax1, t_goal_yr, '--', 'Color',[0.3 0.9 0.3], 'LineWidth',1.5, ...
    'Label',sprintf('  %d-year goal',t_goal_yr), 'LabelHorizontalAlignment','left','FontSize',8);
ax1.YColor = [0.35 0.75 1.00];
ylabel(ax1,'Total time to 100 AU (years)','Color',[0.35 0.75 1.00]);

yyaxis(ax1,'right');
plot(ax1, dv_sweep, vinf_sweep, '-', 'Color',[1.00 0.65 0.25], 'LineWidth',2.0);
ax1.YColor = [1.00 0.65 0.25];
ylabel(ax1,'v_\infty post-Oberth (km/s)','Color',[1.00 0.65 0.25]);

if ~isnan(dv_for_10yr)
    yyaxis(ax1,'left');
    xline(ax1, dv_for_10yr, '--', 'Color',[0.3 0.9 0.3], 'LineWidth',1.2, ...
        'Label',sprintf('  %.1f km/s for 10yr',dv_for_10yr), ...
        'FontSize',8,'LabelVerticalAlignment','bottom');
end

xlabel(ax1,'Oberth \DeltaV at perihelion (km/s)','Color',txtCol);
title(ax1, sprintf('Solar Oberth Trade  |  perihelion = %.4f AU (%.1f R_\\odot)\n v_{peri} = %.1f km/s   v_{esc} = %.1f km/s', ...
    r_peri_actual/AU, r_peri_actual/R_sun, v_peri, v_esc_peri), ...
    'Color',txtCol,'FontSize',9);
noteStr = sprintf('Best trajectory:\n  Launch: %s\n  Flyby:  %s\n  C3 = %.1f km^2/s^2', ...
    datestr(jd_launch_best-1721058.5,'mmm yyyy'), ...
    datestr(jd_jup_best-1721058.5,'mmm yyyy'), C3_grid(iL_best,iJ_best));
text(ax1,0.98,0.98,noteStr,'Units','normalized','VerticalAlignment','top', ...
    'HorizontalAlignment','right','Color',txtCol,'FontSize',8,'FontName','Monospaced', ...
    'BackgroundColor',[0.10 0.10 0.16 0.85],'Margin',4);

%% Plot 2: Full solar system trajectory -----------------------------------
fig2 = figure('Name','100 AU Mission Trajectory','NumberTitle','off', ...
    'Color',bgCol,'Position',[100 100 880 860]);
ax2 = axes('Parent',fig2,'Color',bgCol,'XColor',axCol,'YColor',axCol, ...
    'GridColor',gridC,'Box','on');
hold(ax2,'on');  grid(ax2,'on');  axis(ax2,'equal');

for body = {bodies.Earth, bodies.Jupiter}
    b = body{1};
    plot(ax2, (b.a/AU)*cos(th), (b.a/AU)*sin(th), '-', 'Color',[0.4 0.4 0.5 0.3], 'LineWidth',0.7);
end

npts   = 1500;
r_leg1 = zeros(2, npts);
for k = 1:npts
    try
        [rk,~] = keplerPropagate(r_earth_dep, v1_best, ...
            (k-1)/(npts-1)*tof_leg1_days*86400, muSun);
        r_leg1(:,k) = rk(1:2)/AU;
    catch
        r_leg1(:,k) = NaN;
    end
end
hL1 = plot(ax2, r_leg1(1,:), r_leg1(2,:), '-', 'Color',[0.35 0.75 1.00], 'LineWidth',2.0);
hL2 = plot(ax2, xy_in(1,:),  xy_in(2,:),  '--','Color',[1.00 0.65 0.25], 'LineWidth',2.0);
hL3 = plot(ax2, xy_esc(1,:), xy_esc(2,:), '-', 'Color',[1.00 0.35 0.35], 'LineWidth',2.5);
quiver(ax2, xy_esc(1,end-1), xy_esc(2,end-1), ...
    xy_esc(1,end)-xy_esc(1,end-1), xy_esc(2,end)-xy_esc(2,end-1), ...
    0, 'Color',[1.00 0.35 0.35], 'LineWidth',2.0, 'MaxHeadSize',3.0);

plot(ax2, r_earth_dep(1)/AU, r_earth_dep(2)/AU, 'o','MarkerSize',8, ...
    'MarkerFaceColor',[0.20 0.45 0.75],'MarkerEdgeColor','w');
text(ax2, r_earth_dep(1)/AU+0.1, r_earth_dep(2)/AU, ...
    sprintf('  Earth\n  %s',datestr(jd_launch_best-1721058.5,'mmm yyyy')), ...
    'Color',txtCol,'FontSize',7);
plot(ax2, r_jup_arr(1)/AU, r_jup_arr(2)/AU, 'o','MarkerSize',12, ...
    'MarkerFaceColor',[0.75 0.60 0.45],'MarkerEdgeColor','w');
text(ax2, r_jup_arr(1)/AU, r_jup_arr(2)/AU+0.4, ...
    sprintf('Jupiter flyby\n%s',datestr(jd_jup_best-1721058.5,'mmm yyyy')), ...
    'Color',txtCol,'FontSize',7,'HorizontalAlignment','center');
plot(ax2, r_peri_vec(1)/AU, r_peri_vec(2)/AU, 'y*', ...
    'MarkerSize',11,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor','w');
text(ax2, r_peri_vec(1)/AU+0.12, r_peri_vec(2)/AU, ...
    sprintf(' perihelion\n %.4f AU',r_peri_actual/AU),'Color',txtCol,'FontSize',7);

fill(ax2, r_sun_AU*cos(th), r_sun_AU*sin(th), [1.00 0.85 0.20], ...
    'EdgeColor',[1.00 0.95 0.50],'LineWidth',0.8);
plot(ax2, (r_peri_min/AU)*cos(th), (r_peri_min/AU)*sin(th), '--', ...
    'Color',[1 0.4 0.1 0.6],'LineWidth',1.0);
text(ax2, r_peri_min/AU+0.03, 0, ' PSP limit','Color',[1 0.5 0.2],'FontSize',7);

xlim(ax2,[-6.5 6.5]);  ylim(ax2,[-6.5 6.5]);
xlabel(ax2,'x (AU, ecliptic J2000)','Color',txtCol);
ylabel(ax2,'y (AU, ecliptic J2000)','Color',txtCol);
title(ax2, sprintf('100 AU Mission Trajectory\nLaunch %s  |  Jupiter flyby %s  |  Oberth \\DeltaV = %.1f km/s  |  v_\\infty = %.1f km/s', ...
    datestr(jd_launch_best-1721058.5,'mmm yyyy'), ...
    datestr(jd_jup_best-1721058.5,'mmm yyyy'), dv_oberth, vinf_grid(iL_best,iJ_best)), ...
    'Color',txtCol,'FontSize',10);
legend([hL1 hL2 hL3],{'Earth\rightarrowJupiter','Jupiter\rightarrowPerihelion','Post-Oberth escape'}, ...
    'Location','northwest','TextColor',txtCol,'FontSize',7, ...
    'Color',[0.10 0.10 0.16],'EdgeColor',axCol);

%% Plot 3: Jupiter flyby geometry (perifocal frame) -----------------------
R_JUP     = r_jup_body;
e_jup_fly = 1 + rp_fly * vinf_arr_mag^2 / mu_jup;
p_jup_fly = rp_fly * (1 + e_jup_fly);
nu_asym_j = acos(-1/e_jup_fly);
nu_fly    = linspace(-nu_asym_j*0.90, nu_asym_j*0.90, 1200);
r_fly     = p_jup_fly ./ (1 + e_jup_fly * cos(nu_fly));
xf        = r_fly .* cos(nu_fly) / R_JUP;
yf        = r_fly .* sin(nu_fly) / R_JUP;

fig3 = figure('Name','Jupiter Flyby Geometry','NumberTitle','off', ...
    'Color',bgCol,'Position',[120 120 760 760]);
ax3 = axes('Parent',fig3,'Color',bgCol,'XColor',axCol,'YColor',axCol, ...
    'GridColor',gridC,'Box','on');
hold(ax3,'on');  grid(ax3,'on');  axis(ax3,'equal');

% SOI
r_soi_rj = bodies.Jupiter.soi / R_JUP;
hSOI = plot(ax3, r_soi_rj*cos(th), r_soi_rj*sin(th), '--', ...
    'Color',[0.5 0.5 0.6 0.6],'LineWidth',1.0);
text(ax3, 0, r_soi_rj, '  SOI','Color',[0.6 0.6 0.7],'FontSize',7, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

% Jupiter disk
fill(ax3, cos(th), sin(th), [0.75 0.60 0.45]*0.55, 'EdgeColor','none');
hJUP = plot(ax3, cos(th), sin(th), '-','Color',[0.85 0.70 0.50],'LineWidth',1.5);

% Flyby arc: dashed inbound, solid outbound
imid = floor(numel(nu_fly)/2) + 1;
hFin  = plot(ax3, xf(1:imid),   yf(1:imid),   '--','Color',[1.00 0.65 0.25],'LineWidth',2.0);
hFout = plot(ax3, xf(imid:end), yf(imid:end), '-', 'Color',[1.00 0.35 0.35],'LineWidth',2.0);
quiver(ax3, xf(end-1), yf(end-1), xf(end)-xf(end-1), yf(end)-yf(end-1), ...
    0,'Color',[1.00 0.35 0.35],'LineWidth',1.5,'MaxHeadSize',3.0);

% Periapsis
hPeri3 = plot(ax3, rp_fly/R_JUP, 0, 'y*','MarkerSize',12, ...
    'MarkerFaceColor',[1 1 0],'MarkerEdgeColor','w');
text(ax3, rp_fly/R_JUP*1.05, 0, sprintf('  r_p = %.2f R_J',rp_fly/R_JUP), ...
    'Color',txtCol,'FontSize',8,'VerticalAlignment','middle');

xlim(ax3, [-1 1]*max(abs(xf))*1.12);
ylim(ax3, [-1 1]*max(abs(yf))*1.12);
xlabel(ax3,'Perifocal  x  (R_{Jup},  periapsis along +x)','Color',txtCol);
ylabel(ax3,'Perifocal  y  (R_{Jup})','Color',txtCol);
title(ax3, sprintf('Jupiter Flyby Geometry\nv_\\infty = %.2f km/s  |  turn angle = %.1f°  |  r_p = %.2f R_{Jup}', ...
    vinf_arr_mag, delta_fly*180/pi, rp_fly/R_JUP), 'Color',txtCol,'FontSize',10);
legend([hSOI hJUP hFin hFout hPeri3], ...
    {'SOI','Jupiter','Inbound','Outbound','Periapsis'}, ...
    'Location','northwest','TextColor',txtCol,'FontSize',8, ...
    'Color',[0.10 0.10 0.16],'EdgeColor',axCol);

%% Plot 4: Solar Oberth maneuver close-up ---------------------------------
zoom_au = 0.22;   % half-width of zoom window (AU)

fig4 = figure('Name','Solar Oberth Maneuver','NumberTitle','off', ...
    'Color',bgCol,'Position',[140 140 760 760]);
ax4 = axes('Parent',fig4,'Color',bgCol,'XColor',axCol,'YColor',axCol, ...
    'GridColor',gridC,'Box','on');
hold(ax4,'on');  grid(ax4,'on');  axis(ax4,'equal');

% Sun
fill(ax4, r_sun_AU*cos(th), r_sun_AU*sin(th), [1.00 0.85 0.20], ...
    'EdgeColor',[1.00 0.95 0.50],'LineWidth',0.8);
hSun4 = plot(ax4, r_sun_AU*cos(th), r_sun_AU*sin(th), '-', ...
    'Color',[1.00 0.95 0.50],'LineWidth',1.0);
text(ax4, 0, r_sun_AU*2.5, 'Sun','Color',[1.00 0.85 0.20],'FontSize',9, ...
    'HorizontalAlignment','center');

% PSP thermal limit
hPSP4 = plot(ax4, (r_peri_min/AU)*cos(th), (r_peri_min/AU)*sin(th), '--', ...
    'Color',[1.00 0.45 0.10 0.8],'LineWidth',1.2);

% Inbound arc: show points within zoom window
mask_in = (abs(xy_in(1,:)) <= zoom_au) & (abs(xy_in(2,:)) <= zoom_au);
if any(mask_in)
    seg = find(mask_in);
    hIn4 = plot(ax4, xy_in(1,seg), xy_in(2,seg), '--', ...
        'Color',[1.00 0.65 0.25],'LineWidth',2.0);
else
    hIn4 = plot(ax4, NaN, NaN,'--','Color',[1.00 0.65 0.25],'LineWidth',2.0);
end

% Post-Oberth escape arc: show points within zoom window
mask_esc = (abs(xy_esc(1,:)) <= zoom_au) & (abs(xy_esc(2,:)) <= zoom_au);
if any(mask_esc)
    seg_e = find(mask_esc);
    hEsc4 = plot(ax4, xy_esc(1,seg_e), xy_esc(2,seg_e), '-', ...
        'Color',[1.00 0.35 0.35],'LineWidth',2.5);
    if numel(seg_e) > 1
        last2 = seg_e(end-1:end);
        quiver(ax4, xy_esc(1,last2(1)), xy_esc(2,last2(1)), ...
            xy_esc(1,last2(2))-xy_esc(1,last2(1)), ...
            xy_esc(2,last2(2))-xy_esc(2,last2(1)), ...
            0,'Color',[1.00 0.35 0.35],'LineWidth',2.0,'MaxHeadSize',3.0);
    end
else
    hEsc4 = plot(ax4, NaN, NaN,'-','Color',[1.00 0.35 0.35],'LineWidth',2.5);
end

% Perihelion marker
hP4 = plot(ax4, r_peri_vec(1)/AU, r_peri_vec(2)/AU, 'y*', ...
    'MarkerSize',14,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor','w');
text(ax4, r_peri_vec(1)/AU, r_peri_vec(2)/AU - zoom_au*0.06, ...
    sprintf('Perihelion\n%.4f AU  (%.1f R_{\\odot})', r_peri_actual/AU, r_peri_actual/R_sun), ...
    'Color',txtCol,'FontSize',7.5,'HorizontalAlignment','center','VerticalAlignment','top');

% ΔV arrow (tangential at perihelion)
dv_scale = zoom_au * 0.14;
hDV4 = quiver(ax4, r_peri_vec(1)/AU, r_peri_vec(2)/AU, ...
    v_peri_dir(1)*dv_scale, v_peri_dir(2)*dv_scale, ...
    0,'Color',[0.30 1.00 0.30],'LineWidth',2.5,'MaxHeadSize',2.5);
text(ax4, r_peri_vec(1)/AU + v_peri_dir(1)*dv_scale*1.25, ...
         r_peri_vec(2)/AU + v_peri_dir(2)*dv_scale*1.25, ...
    sprintf(' \\DeltaV = %.1f km/s',dv_oberth), ...
    'Color',[0.30 1.00 0.30],'FontSize',8);

xlim(ax4,[-zoom_au zoom_au]);  ylim(ax4,[-zoom_au zoom_au]);
xlabel(ax4,'x (AU, ecliptic J2000)','Color',txtCol);
ylabel(ax4,'y (AU, ecliptic J2000)','Color',txtCol);
title(ax4, sprintf('Solar Oberth Maneuver\nr_{peri} = %.4f AU (%.1f R_\\odot)  |  v_{peri} = %.1f km/s  |  v_{esc} = %.1f km/s  |  v_\\infty = %.1f km/s', ...
    r_peri_actual/AU, r_peri_actual/R_sun, v_peri, v_esc_peri, vinf_actual), ...
    'Color',txtCol,'FontSize',9.5);
legend([hSun4 hPSP4 hIn4 hEsc4 hP4 hDV4], ...
    {'Sun','PSP thermal limit','Inbound arc','Post-Oberth escape','Perihelion','\DeltaV'}, ...
    'Location','northwest','TextColor',txtCol,'FontSize',8, ...
    'Color',[0.10 0.10 0.16],'EdgeColor',axCol);

fprintf('Done.\n');

%% ---- Local helpers -----------------------------------------------------
function [r1, v1] = keplerPropagate(r0, v0, dt, mu)
    r0m    = norm(r0);
    v0m    = norm(v0);
    alpha  = 2/r0m - v0m^2/mu;
    sigma0 = dot(r0, v0) / sqrt(mu);
    chi    = sqrt(mu) * abs(alpha) * dt;
    for k = 1:50
        psi       = chi^2 * alpha;
        [C, S]    = stumpffCS(psi);
        r_now     = chi^2*C + sigma0*chi*(1-psi*S) + r0m*(1-psi*C);
        dchi      = (sqrt(mu)*dt - chi^3*S - sigma0*chi^2*C - r0m*chi*(1-psi*S)) / r_now;
        chi       = chi + dchi;
        if abs(dchi) < 1e-10, break; end
    end
    psi  = chi^2 * alpha;
    [C, S] = stumpffCS(psi);
    f  =  1 - chi^2*C / r0m;
    g  =  dt - chi^3*S / sqrt(mu);
    r1 =  f*r0 + g*v0;
    r1m = norm(r1);
    df =  sqrt(mu)*chi*(psi*S - 1) / (r0m*r1m);
    dg =  1 - chi^2*C / r1m;
    v1 =  df*r0 + dg*v0;
end

function [C, S] = stumpffCS(psi)
    if psi > 1e-6
        sp = sqrt(psi);
        C  = (1 - cos(sp))  / psi;
        S  = (sp - sin(sp)) / (psi*sp);
    elseif psi < -1e-6
        sp = sqrt(-psi);
        C  = (1 - cosh(sp)) / psi;
        S  = (sinh(sp) - sp) / ((-psi)*sp);
    else
        C  = 0.5;
        S  = 1/6;
    end
end
