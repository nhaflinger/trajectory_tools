function result = gravityAssist(body, v_inf_in_vec, v_inf_out_vec, options)
%GRAVITYASSIST  Analyse a gravity-assist flyby using patched-conic v∞ vectors.
%
%   result = gravityAssist(body, v_inf_in_vec, v_inf_out_vec)
%   result = gravityAssist(body, v_inf_in_vec, v_inf_out_vec, options)
%
%   Given the incoming and outgoing hyperbolic-excess velocity vectors in
%   the heliocentric frame (spacecraft velocity minus planet velocity at
%   the flyby body), computes the required flyby geometry and feasibility.
%
%   A "free" gravity assist requires |v∞_in| = |v∞_out|.  When Lambert
%   solutions for adjacent legs give different magnitudes the flyby is
%   "powered": a burn at periapsis bridges the speed difference via the
%   Oberth effect.
%
%   Inputs:
%     body              - body struct from constants() (.mu, .radius, .name)
%     v_inf_in_vec      - 3×1 incoming v∞ vector (km/s, heliocentric frame)
%     v_inf_out_vec     - 3×1 outgoing v∞ vector (km/s, heliocentric frame)
%     options
%       .atmosphereAltitude  altitude added to body radius as minimum safe
%                            periapsis (km) [default 0]
%       .flybyAltitude       periapsis altitude to use when computing powered-
%                            flyby ΔV (km) [default 300]
%
%   Output fields:
%     .v_inf_in      incoming hyperbolic excess speed (km/s)
%     .v_inf_out     outgoing hyperbolic excess speed (km/s)
%     .deflection    turning angle δ between v∞_in and v∞_out (deg)
%     .r_periapsis   closest-approach radius required for deflection (km)
%                    (= body.radius + flybyAltitude for powered flybys)
%     .altitude      closest-approach altitude above body surface (km)
%     .isFeasible    true when periapsis ≥ body.radius + atmosphereAltitude
%     .isPowered     true when |v∞_in| ≠ |v∞_out|
%     .dvPowered     ΔV at periapsis for powered component (km/s); 0 for free GA
%     .maxDeflection maximum achievable deflection at minimum safe periapsis (deg)
%     .v_inf_in_vec  3×1 incoming v∞ vector (stored for plotting)
%     .v_inf_out_vec 3×1 outgoing v∞ vector

if nargin < 4, options = struct(); end
if ~isfield(options, 'atmosphereAltitude'), options.atmosphereAltitude = 0;   end
if ~isfield(options, 'flybyAltitude'),      options.flybyAltitude      = 300; end

v_inf_in_vec  = v_inf_in_vec(:);
v_inf_out_vec = v_inf_out_vec(:);

v_in  = norm(v_inf_in_vec);
v_out = norm(v_inf_out_vec);
mu    = body.mu;
r_min = body.radius + options.atmosphereAltitude;

% --- Deflection angle ---
cos_delta = dot(v_inf_in_vec, v_inf_out_vec) / (v_in * v_out);
cos_delta = max(-1.0, min(1.0, cos_delta));
delta_rad = acos(cos_delta);
delta_deg = rad2deg(delta_rad);

% --- Classify flyby type ---
isPowered = abs(v_in - v_out) > 1e-4;   % 0.1 m/s threshold

if ~isPowered
    % Pure gravity assist: sin(δ/2) = 1 / (1 + r_p·v∞²/μ)
    %   → r_p = (μ/v∞²) · (1/sin(δ/2) − 1)
    v_inf    = 0.5 * (v_in + v_out);
    sin_half = sin(delta_rad / 2);
    if sin_half < 1e-10
        r_periapsis = Inf;
    else
        r_periapsis = (mu / v_inf^2) * (1/sin_half - 1);
    end
    altitude  = r_periapsis - body.radius;
    dvPowered = 0;
else
    % Powered flyby: ΔV applied at specified periapsis altitude
    r_periapsis = body.radius + options.flybyAltitude;
    altitude    = options.flybyAltitude;
    v_p_in  = sqrt(v_in^2  + 2*mu/r_periapsis);
    v_p_out = sqrt(v_out^2 + 2*mu/r_periapsis);
    dvPowered = abs(v_p_out - v_p_in);
end

isFeasible = isfinite(r_periapsis) && (r_periapsis >= r_min);

% --- Maximum achievable deflection at minimum safe periapsis ---
v_ref    = 0.5 * (v_in + v_out);
sin_max  = 1 / (1 + r_min * v_ref^2 / mu);
dmax_deg = rad2deg(2 * asin(min(1.0, sin_max)));

result = struct( ...
    'v_inf_in',      v_in,       ...
    'v_inf_out',     v_out,      ...
    'deflection',    delta_deg,  ...
    'r_periapsis',   r_periapsis, ...
    'altitude',      altitude,   ...
    'isFeasible',    isFeasible, ...
    'isPowered',     isPowered,  ...
    'dvPowered',     dvPowered,  ...
    'maxDeflection', dmax_deg,   ...
    'v_inf_in_vec',  v_inf_in_vec, ...
    'v_inf_out_vec', v_inf_out_vec);
end
