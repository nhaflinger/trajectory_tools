function res = planeChange(orb, delta_i_deg, varargin)
%PLANECHANGE  Plane change maneuver analysis for a circular Earth orbit.
%
%   res = planeChange(orb, delta_i_deg)
%   res = planeChange(orb, delta_i_deg, 'Type', 'combined', 'AltFinal_km', alt2)
%
%   orb          - earthOrbit struct (circular orbit assumed; uses orb.a)
%   delta_i_deg  - magnitude of inclination change in degrees (>= 0)
%
%   Options:
%     'Type'         - 'pure' (default) or 'combined'
%     'AltFinal_km'  - final altitude in km (required for 'combined')
%
%   'pure': pure plane change only, no altitude change.
%     dv = 2 * v_circ * sin(delta_i/2)
%
%   'combined': simultaneous plane change + Hohmann altitude change at apogee
%     of the Hohmann transfer orbit.
%     dv_combined = sqrt(v_t_apo^2 + v_circ2^2 - 2*v_t_apo*v_circ2*cos(delta_i))
%
%   Output struct fields:
%     dv_km_s          - DV for the maneuver
%     dv_separate_km_s - DV if done separately (combined mode only: hohmann dv2 + pure)
%     savings_km_s     - dv_separate - dv_combined  (combined mode only)
%     delta_i_deg
%     type             - 'pure' or 'combined'
%     alt_km           - altitude (pure) or {alt1_km, alt2_km} style fields (combined)

mu_E = 398600.4418;   % km^3/s^2
R_E  = 6378.1363;     % km

%% ── Parse options ────────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'Type',        'pure', @(x) any(strcmpi(x, {'pure','combined'})));
addParameter(p, 'AltFinal_km', [],     @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});
type_str   = lower(p.Results.Type);
alt2       = p.Results.AltFinal_km;

%% ── Extract orbit parameters ─────────────────────────────────────────────────
r1     = orb.a;           % circular orbit radius (km)
alt1   = r1 - R_E;
v_circ = sqrt(mu_E / r1); % circular orbital speed

%% ── Compute maneuver ─────────────────────────────────────────────────────────
switch type_str

    case 'pure'
        dv = 2 * v_circ * sind(delta_i_deg / 2);

        res = struct( ...
            'type',             'pure',     ...
            'delta_i_deg',      delta_i_deg, ...
            'alt_km',           alt1,        ...
            'dv_km_s',          dv);

        fprintf('\n=== Pure Plane Change ===\n');
        fprintf('  Altitude         : %.2f km\n', alt1);
        fprintf('  v_circular       : %.4f km/s\n', v_circ);
        fprintf('  Delta-i          : %.4f deg\n', delta_i_deg);
        fprintf('  DV               : %.4f km/s  (%.2f m/s)\n', dv, dv*1000);

    case 'combined'
        if isempty(alt2)
            error('planeChange: ''AltFinal_km'' must be supplied for combined type.');
        end
        r2       = R_E + alt2;
        a_t      = (r1 + r2) / 2;

        % Transfer orbit velocity at the apogee (far end from burn point)
        % Burn point is at r1 (lower orbit), transfer goes to r2 (higher orbit),
        % plane change applied at r2 (apogee of Hohmann transfer).
        v_t_apo  = sqrt(mu_E * (2/r2 - 1/a_t));   % transfer velocity at r2
        v_circ2  = sqrt(mu_E / r2);                 % circular velocity at r2

        dv_combined = sqrt(v_t_apo^2 + v_circ2^2 - ...
                           2*v_t_apo*v_circ2*cosd(delta_i_deg));

        % Separate comparison: Hohmann dv2 (circularisation) + pure plane change at r2
        dv_hohmann_dv2 = abs(v_circ2 - v_t_apo);
        dv_pure_r2     = 2 * v_circ2 * sind(delta_i_deg / 2);
        dv_separate    = dv_hohmann_dv2 + dv_pure_r2;

        savings = dv_separate - dv_combined;

        res = struct( ...
            'type',              'combined',    ...
            'delta_i_deg',       delta_i_deg,   ...
            'alt1_km',           alt1,           ...
            'alt2_km',           alt2,           ...
            'dv_km_s',           dv_combined,    ...
            'dv_separate_km_s',  dv_separate,    ...
            'savings_km_s',      savings);

        fprintf('\n=== Combined Plane Change + Altitude Change ===\n');
        fprintf('  Initial altitude  : %.2f km  (v = %.4f km/s)\n', alt1, v_circ);
        fprintf('  Final altitude    : %.2f km  (v = %.4f km/s)\n', alt2, v_circ2);
        fprintf('  Delta-i           : %.4f deg\n', delta_i_deg);
        fprintf('  Transfer apo vel  : %.4f km/s\n', v_t_apo);
        fprintf('  DV combined       : %.4f km/s  (%.2f m/s)\n', dv_combined, dv_combined*1000);
        fprintf('  DV separate       : %.4f km/s  (%.2f m/s)\n', dv_separate, dv_separate*1000);
        if savings >= 0
            fprintf('  Savings (combined): %.4f km/s  (%.2f m/s)  -- combined is BETTER\n', savings, savings*1000);
        else
            fprintf('  Savings (combined): %.4f km/s  (%.2f m/s)  -- separate is BETTER\n', savings, savings*1000);
        end

    otherwise
        error('planeChange: unknown Type ''%s''. Use ''pure'' or ''combined''.', type_str);
end

end
