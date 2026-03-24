function res = deorbitBurn(orb_or_alt, varargin)
%DEORBITBURN  Retrograde deorbit burn analysis from a circular Earth orbit.
%
%   res = deorbitBurn(orb)
%   res = deorbitBurn(orb, 'TargetPeriAlt_km', 80)
%   res = deorbitBurn(alt_km)
%   res = deorbitBurn(alt_km, 'TargetPeriAlt_km', 80)
%
%   orb_or_alt         - earthOrbit struct OR scalar circular altitude (km)
%   'TargetPeriAlt_km' - desired perigee altitude after burn (km). Default: 80 km
%                        (standard reentry interface)
%
%   Physics: lowers perigee from the circular orbit radius to target_alt.
%     a_deorbit = (r_circ + r_peri) / 2
%     dv = v_circ - v_transfer_at_apogee  (retrograde, positive value)
%
%   Output struct fields:
%     dv_km_s              - deorbit delta-V (km/s)
%     dv_m_s               - deorbit delta-V (m/s)
%     alt_km               - initial circular orbit altitude (km)
%     target_peri_alt_km   - target perigee altitude (km)
%     tof_to_peri_s        - time of flight from burn to perigee (half-period, s)
%     tof_to_peri_min      - same in minutes
%     a_deorbit_km         - semi-major axis of deorbit ellipse (km)
%     compliant_25yr       - logical; true if alt_km <= 600 km (passive reentry rule)

mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km

%% ── Parse input ──────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'TargetPeriAlt_km', 80, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});
target_peri_alt = p.Results.TargetPeriAlt_km;

if isstruct(orb_or_alt)
    r_circ = orb_or_alt.a;        % circular orbit radius (km)
    alt_km = r_circ - R_E;
elseif isnumeric(orb_or_alt) && isscalar(orb_or_alt)
    alt_km = orb_or_alt;
    r_circ = R_E + alt_km;
else
    error('deorbitBurn: first argument must be an earthOrbit struct or scalar altitude in km.');
end

%% ── Validate ─────────────────────────────────────────────────────────────────
if target_peri_alt >= alt_km
    error('deorbitBurn: target perigee altitude (%.1f km) must be below the orbit altitude (%.1f km).', ...
          target_peri_alt, alt_km);
end
if target_peri_alt < 0
    error('deorbitBurn: target perigee altitude cannot be negative (%.1f km).', target_peri_alt);
end

%% ── Physics ──────────────────────────────────────────────────────────────────
r_peri     = R_E + target_peri_alt;
a_deorbit  = (r_circ + r_peri) / 2;

v_circ         = sqrt(mu_E / r_circ);
v_trans_at_apo = sqrt(mu_E * (2/r_circ - 1/a_deorbit));

dv = v_circ - v_trans_at_apo;   % positive retrograde burn

tof_to_peri = pi * sqrt(a_deorbit^3 / mu_E);   % half-period

compliant_25yr = (alt_km <= 600);

%% ── Output ───────────────────────────────────────────────────────────────────
res = struct( ...
    'dv_km_s',            dv,                ...
    'dv_m_s',             dv * 1000,         ...
    'alt_km',             alt_km,            ...
    'target_peri_alt_km', target_peri_alt,   ...
    'tof_to_peri_s',      tof_to_peri,       ...
    'tof_to_peri_min',    tof_to_peri / 60,  ...
    'a_deorbit_km',       a_deorbit,         ...
    'compliant_25yr',     compliant_25yr);

fprintf('\n=== Deorbit Burn Analysis ===\n');
fprintf('  Circular orbit alt  : %.2f km  (r = %.2f km)\n', alt_km, r_circ);
fprintf('  Target perigee alt  : %.2f km  (r = %.2f km)\n', target_peri_alt, r_peri);
fprintf('  Deorbit ellipse SMA : %.2f km\n', a_deorbit);
fprintf('  Circular speed      : %.4f km/s\n', v_circ);
fprintf('  Transfer speed (apo): %.4f km/s\n', v_trans_at_apo);
fprintf('  DV (retrograde)     : %.4f km/s  (%.2f m/s)\n', dv, dv*1000);
fprintf('  TOF to perigee      : %.2f s  (%.2f min)\n', tof_to_peri, tof_to_peri/60);
if compliant_25yr
    fprintf('  25-year rule        : COMPLIANT (alt <= 600 km)\n');
else
    fprintf('  25-year rule        : NON-COMPLIANT (alt > 600 km, active deorbit may be needed)\n');
end

end
