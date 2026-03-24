function result = cwRendezvous(dr0, dv0, n, tf, varargin)
%CWRENDEZVOUS  Two-impulse CW rendezvous using the Clohessy-Wiltshire STM.
%
%   result = cwRendezvous(dr0, dv0, n, tf)
%   result = cwRendezvous(dr0, dv0, n, tf, drf, dvf)
%
%   Computes two impulsive delta-V manoeuvres to transfer a deputy from an
%   initial relative state (dr0, dv0) to a desired final relative state
%   (drf, dvf) in transfer time tf, using the Clohessy-Wiltshire equations.
%
%   VELOCITY CONVENTION: dv0, dvf, and output dv1/dv2 are all CW rotating-frame
%   derivatives (d/dt in the rotating LVLH frame). See cwPropagate for details.
%
%   LVLH frame: R = radial out, S = along-track, W = orbit normal.
%
%   SINGULARITY WARNING: Phi_rv is singular when tf = k*(T/2) for any integer k
%   (i.e. at half-period multiples). The W component becomes uncontrollable at
%   tf = k*T, and the full matrix is near-singular at tf = k*T/2. Choose tf
%   away from these values — e.g. 0.4*T, 0.6*T, 0.75*T are all well-conditioned.
%
%   Algorithm:
%     The CW STM relates initial and final states via:
%       drf = Phi_rr * dr0 + Phi_rv * (dv0 + Dv1)
%     Solving for Dv1:
%       Dv1 = Phi_rv \ (drf - Phi_rr * dr0) - dv0
%     After free drift from dr0 with (dv0 + Dv1), the velocity at tf is:
%       dv_arrive = Phi_vr * dr0 + Phi_vv * (dv0 + Dv1)
%     Second burn nulls the residual:
%       Dv2 = dvf - dv_arrive
%
%   Inputs:
%     dr0  - initial relative position (km), 3x1, LVLH [R; S; W]
%     dv0  - initial relative velocity (km/s), 3x1, CW rotating-frame
%     n    - chief mean motion (rad/s)
%     tf   - transfer time (s)
%     drf  - desired final relative position (km), 3x1 (default: [0;0;0])
%     dvf  - desired final relative velocity (km/s), 3x1 (default: [0;0;0])
%
%   Output struct fields:
%     dv1          - first burn delta-V vector (km/s), 3x1
%     dv2          - second burn delta-V vector (km/s), 3x1
%     dv1_mag      - first burn magnitude (m/s)
%     dv2_mag      - second burn magnitude (m/s)
%     total_dv_m_s - total delta-V (m/s)
%     tf           - transfer time (s)
%     r_hist       - relative position trajectory 3x100 (km)
%     v_hist       - relative velocity trajectory 3x100 (km/s)
%     t_hist       - time vector 1x100 (s)

% Parse optional arguments
if numel(varargin) >= 1
    drf = varargin{1}(:);
else
    drf = [0; 0; 0];
end
if numel(varargin) >= 2
    dvf = varargin{2}(:);
else
    dvf = [0; 0; 0];
end

dr0 = dr0(:);
dv0 = dv0(:);

tau = n * tf;
s   = sin(tau);
c   = cos(tau);

% ── CW STM at tf ────────────────────────────────────────────────────────────
Phi_rr = [4 - 3*c,          0,  0;
           6*(s - tau),     1,  0;
           0,                0,  c];

Phi_rv = [ s/n,              2*(1 - c)/n,    0;
           2*(c - 1)/n,     (4*s - 3*tau)/n, 0;
           0,                0,               s/n];

Phi_vr = [ 3*n*s,           0,  0;
           6*n*(c - 1),     0,  0;
           0,                0,  -n*s];

Phi_vv = [ c,       2*s,       0;
           -2*s,    4*c - 3,   0;
           0,        0,         c];

% ── Singularity check ────────────────────────────────────────────────────────
cond_rv = cond(Phi_rv);
if cond_rv > 1e8
    warning('cwRendezvous: Phi_rv is near-singular (cond = %.3e). Transfer time tf = %.1f s is likely near a resonant multiple of the orbital period (T = %.1f s). Results may be inaccurate.', ...
        cond_rv, tf, 2*pi/n);
end

% ── Solve for first burn ─────────────────────────────────────────────────────
% drf = Phi_rr * dr0 + Phi_rv * (dv0 + dv1)
% Phi_rv * dv1 = drf - Phi_rr * dr0 - Phi_rv * dv0
rhs_1 = drf - Phi_rr * dr0 - Phi_rv * dv0;
dv1   = Phi_rv \ rhs_1;

% Velocity after first burn
dv_post_burn1 = dv0 + dv1;

% ── Velocity at arrival (before second burn) ─────────────────────────────────
dv_arrive = Phi_vr * dr0 + Phi_vv * dv_post_burn1;

% ── Second burn ──────────────────────────────────────────────────────────────
dv2 = dvf - dv_arrive;

% ── Magnitudes ───────────────────────────────────────────────────────────────
dv1_mag      = norm(dv1) * 1e3;      % m/s
dv2_mag      = norm(dv2) * 1e3;      % m/s
total_dv_m_s = dv1_mag + dv2_mag;

% ── Trajectory between burns ─────────────────────────────────────────────────
N_traj = 100;
t_hist = linspace(0, tf, N_traj);
[r_hist, v_hist] = cwPropagate(dr0, dv_post_burn1, n, t_hist);

% ── Build output struct ───────────────────────────────────────────────────────
result = struct( ...
    'dv1',          dv1,          ...
    'dv2',          dv2,          ...
    'dv1_mag',      dv1_mag,      ...
    'dv2_mag',      dv2_mag,      ...
    'total_dv_m_s', total_dv_m_s, ...
    'tf',           tf,           ...
    'r_hist',       r_hist,       ...
    'v_hist',       v_hist,       ...
    't_hist',       t_hist);

% ── Print summary ─────────────────────────────────────────────────────────────
T_chief = 2*pi / n;
fprintf('=== cwRendezvous Summary ===\n');
fprintf('  Transfer time : %.1f s (%.4f orbits)\n', tf, tf / T_chief);
fprintf('  Initial pos   : [%.3f, %.3f, %.3f] km  (R,S,W)\n', dr0(1), dr0(2), dr0(3));
fprintf('  Target  pos   : [%.3f, %.3f, %.3f] km  (R,S,W)\n', drf(1), drf(2), drf(3));
fprintf('  dv1 vector    : [%.4f, %.4f, %.4f] km/s\n', dv1(1), dv1(2), dv1(3));
fprintf('  dv2 vector    : [%.4f, %.4f, %.4f] km/s\n', dv2(1), dv2(2), dv2(3));
fprintf('  |dv1|         : %.2f m/s\n', dv1_mag);
fprintf('  |dv2|         : %.2f m/s\n', dv2_mag);
fprintf('  Total dV      : %.2f m/s\n', total_dv_m_s);
fprintf('============================\n');

end
