function [r_hist, v_hist] = cwPropagate(r0, v0, n, t)
%CWPROPAGATE  Analytical Clohessy-Wiltshire (CW) state propagation.
%
%   [r_hist, v_hist] = cwPropagate(r0, v0, n, t)
%
%   Propagates a deputy spacecraft's relative state in the LVLH (Local
%   Vertical / Local Horizontal) rotating frame using the Clohessy-Wiltshire
%   state transition matrix (STM). Valid for near-circular chief orbits and
%   small relative separations.
%
%   IMPORTANT — Velocity convention:
%     v0 and v_hist are CW rotating-frame derivatives (d/dt in the rotating
%     frame), NOT inertial relative velocities. The relationship between the
%     two is:
%       v_inertial_relative = v_CW + omega x r_LVLH
%     where omega = n * W_hat = [0; 0; n] in LVLH coordinates.
%     To initialise from an ECI state:
%       v_CW = R_eci2lvlh * (v_deputy_ECI - v_chief_ECI) - cross([0;0;n], r_LVLH)
%
%   LVLH frame convention:
%     R — radial (outward positive)
%     S — along-track (velocity direction positive)
%     W — cross-track (orbit normal, right-hand rule)
%
%   Inputs:
%     r0  - initial relative position, [R; S; W] km, 3x1
%     v0  - initial relative velocity (CW rotating-frame), [dR; dS; dW] km/s, 3x1
%     n   - chief mean motion (rad/s)
%     t   - time vector (s), 1xN or Nx1
%
%   Outputs:
%     r_hist - relative position history, 3xN (km)
%     v_hist - relative velocity history, 3xN (km/s), rotating-frame derivatives

t    = t(:)';    % ensure row vector
N    = numel(t);
r0   = r0(:);
v0   = v0(:);

r_hist = zeros(3, N);
v_hist = zeros(3, N);

for k = 1:N
    tau = n * t(k);    % dimensionless time

    s = sin(tau);
    c = cos(tau);

    % ── State Transition Matrix sub-blocks ──────────────────────────────────
    % Φ_rr: maps r0 → r(t)
    Phi_rr = [4 - 3*c,          0,  0;
               6*(s - tau),     1,  0;
               0,                0,  c];

    % Φ_rv: maps v0 → r(t)  [units: s, so r = Φ_rv * v0 gives km when v0 in km/s]
    Phi_rv = [ s/n,              2*(1 - c)/n,    0;
               2*(c - 1)/n,     (4*s - 3*tau)/n, 0;
               0,                0,               s/n];

    % Φ_vr: maps r0 → v(t)  [units: 1/s, so v = Φ_vr * r0 gives km/s when r0 in km]
    Phi_vr = [ 3*n*s,           0,  0;
               6*n*(c - 1),     0,  0;
               0,                0,  -n*s];

    % Φ_vv: maps v0 → v(t)  [dimensionless]
    Phi_vv = [ c,       2*s,       0;
               -2*s,    4*c - 3,   0;
               0,        0,         c];

    % ── Apply STM ───────────────────────────────────────────────────────────
    r_hist(:, k) = Phi_rr * r0 + Phi_rv * v0;
    v_hist(:, k) = Phi_vr * r0 + Phi_vv * v0;
end

end
