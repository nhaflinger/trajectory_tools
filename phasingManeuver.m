function res = phasingManeuver(orb, delta_phase_deg, n_phasing_orbits)
%PHASINGMANEUVER  Compute a phasing maneuver for a spacecraft in a circular orbit.
%
%   res = phasingManeuver(orb, delta_phase_deg, n_phasing_orbits)
%
%   orb               - earthOrbit struct (circular orbit assumed; uses orb.a, orb.period)
%   delta_phase_deg   - desired phase change in degrees.
%                       Positive = catch up (target is AHEAD of chaser).
%                       Negative = open gap (target is BEHIND chaser).
%   n_phasing_orbits  - number of phasing orbits to complete the maneuver.
%
%   The chaser temporarily enters a phasing orbit whose period differs from
%   the nominal orbit such that after n_phasing_orbits the desired phase offset
%   is achieved.  Two equal-magnitude burns are required (entry and exit).
%
%   Output struct fields:
%     dv_each_km_s         - DV for each burn (km/s)
%     dv_total_km_s        - total DV (2 burns, km/s)
%     a_phase_km           - semi-major axis of phasing orbit (km)
%     alt_peri_phase_km    - perigee altitude of phasing orbit (km)
%     alt_apo_phase_km     - apogee altitude of phasing orbit (km)
%     period_phase_s       - period of phasing orbit (s)
%     period_phase_min     - period of phasing orbit (min)
%     time_to_complete_s   - total wall-clock time for phasing (s)
%     time_to_complete_hr  - total wall-clock time for phasing (hr)
%     n_phasing_orbits     - (echoed)
%     delta_phase_deg      - (echoed)
%     direction            - 'catch_up' or 'open_gap'

mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km

%% ── Extract orbit parameters ─────────────────────────────────────────────────
a       = orb.a;           % nominal orbit SMA (km)
T_nom   = orb.period;      % nominal period (s)
v_circ  = sqrt(mu_E / a);  % circular orbital speed

%% ── Phasing orbit period ─────────────────────────────────────────────────────
% After n_phasing_orbits on the phasing orbit, the satellite has gained
% delta_phase_deg relative to what it would have done on the nominal orbit.
% T_phase * n = T_nom * n - (delta_phase_deg/360) * T_nom
% => T_phase = T_nom * (1 - delta_phase_deg / (360 * n_phasing_orbits))
% Positive delta => shorter period (higher mean motion) => inner orbit => catch up.
T_phase = T_nom * n_phasing_orbits / (n_phasing_orbits - delta_phase_deg/360);
a_phase = (mu_E * (T_phase / (2*pi))^2)^(1/3);

%% ── Check for re-entry ───────────────────────────────────────────────────────
% For a Hohmann-like phasing orbit entered from the current circular orbit,
% the orbit passes through r = a (the nominal orbit radius) at two points,
% so the other extreme is 2*a_phase - a.
r_other = 2*a_phase - a;
alt_other = r_other - R_E;

if r_other < R_E
    error('phasingManeuver: phasing orbit periapsis (%.1f km alt) is below Earth surface. Reduce |delta_phase_deg| or increase n_phasing_orbits.', ...
          alt_other);
end

%% ── DV computation ───────────────────────────────────────────────────────────
v_phase_at_junction = sqrt(mu_E * (2/a - 1/a_phase));
dv_each  = abs(v_phase_at_junction - v_circ);
dv_total = 2 * dv_each;

%% ── Geometry of phasing orbit ────────────────────────────────────────────────
% The phasing orbit shares one apse with the circular orbit radius a.
% For T_phase < T_nom (catch up): phasing orbit is smaller => perigee at r_other, apogee at a
% For T_phase > T_nom (open gap): phasing orbit is larger => perigee at a, apogee at r_other
if a_phase < a
    alt_peri_phase = alt_other;
    alt_apo_phase  = a - R_E;
else
    alt_peri_phase = a - R_E;
    alt_apo_phase  = alt_other;
end

time_total = T_phase * n_phasing_orbits;

if delta_phase_deg > 0
    direction = 'catch_up';
else
    direction = 'open_gap';
end

res = struct( ...
    'dv_each_km_s',        dv_each,          ...
    'dv_total_km_s',       dv_total,         ...
    'a_phase_km',          a_phase,          ...
    'alt_peri_phase_km',   alt_peri_phase,   ...
    'alt_apo_phase_km',    alt_apo_phase,    ...
    'period_phase_s',      T_phase,          ...
    'period_phase_min',    T_phase/60,       ...
    'time_to_complete_s',  time_total,       ...
    'time_to_complete_hr', time_total/3600,  ...
    'n_phasing_orbits',    n_phasing_orbits, ...
    'delta_phase_deg',     delta_phase_deg,  ...
    'direction',           direction);

fprintf('\n=== Phasing Maneuver ===\n');
fprintf('  Nominal orbit SMA   : %.2f km  (alt = %.2f km)\n', a, a - R_E);
fprintf('  Nominal period      : %.2f min\n', T_nom/60);
fprintf('  Phase change        : %+.2f deg  (%s)\n', delta_phase_deg, direction);
fprintf('  Phasing orbits      : %d\n', n_phasing_orbits);
fprintf('  Phasing orbit SMA   : %.2f km\n', a_phase);
fprintf('  Phasing orbit peri  : %.2f km alt\n', alt_peri_phase);
fprintf('  Phasing orbit apo   : %.2f km alt\n', alt_apo_phase);
fprintf('  Phasing period      : %.2f min\n', T_phase/60);
fprintf('  DV (each burn)      : %.4f km/s  (%.2f m/s)\n', dv_each, dv_each*1000);
fprintf('  DV total            : %.4f km/s  (%.2f m/s)\n', dv_total, dv_total*1000);
fprintf('  Time to complete    : %.2f hr  (%.2f min)\n', time_total/3600, time_total/60);

end
