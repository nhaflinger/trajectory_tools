function [v1, v2] = lambertSolver(r1_vec, r2_vec, tof, mu, isLongWay)
%LAMBERTSOLVER  Universal-variable Lambert solver (Bate-Mueller-White §5.3).
%
%   [v1, v2] = lambertSolver(r1_vec, r2_vec, tof, mu)
%   [v1, v2] = lambertSolver(r1_vec, r2_vec, tof, mu, isLongWay)
%
%   Inputs:
%     r1_vec    - 3x1 departure position vector (km)
%     r2_vec    - 3x1 arrival position vector (km)
%     tof       - time of flight (s)
%     mu        - gravitational parameter of central body (km^3/s^2)
%     isLongWay - true for >180 deg transfer (default false)
%
%   Outputs:
%     v1 - 3x1 departure velocity (km/s)
%     v2 - 3x1 arrival velocity (km/s)
%
%   Method — universal variable z parameterises the transfer orbit:
%     A      = dm * sqrt(R1*R2*(1 + cos(dnu)))   dm=+1 short, -1 long
%     y(z)   = R1 + R2 + A*(z*S(z) - 1)/sqrt(C(z))
%     x(z)   = sqrt(y/C(z))
%     TOF(z) = (x^3*S(z) + A*sqrt(y)) / sqrt(mu)
%
%   Velocities from Lagrange f, g, gdot coefficients:
%     f    = 1 - y/R1
%     g    = A*sqrt(y/mu)
%     gdot = 1 - y/R2
%     v1   = (r2 - f*r1) / g
%     v2   = (gdot*r2 - r1) / g
%
%   TOF(z) is monotonically decreasing on (-inf, (2*pi)^2), so a
%   coarse scan locates the bracket and bisection refines it.

if nargin < 5 || isempty(isLongWay)
    isLongWay = false;
end

r1_vec = r1_vec(:);
r2_vec = r2_vec(:);
R1 = norm(r1_vec);
R2 = norm(r2_vec);

% ---- geometry parameter A ----
cos_dnu = dot(r1_vec, r2_vec) / (R1 * R2);
cos_dnu = max(-1.0, min(1.0, cos_dnu));   % numerical clamp
dm = 1 - 2*double(isLongWay);             % +1 short-way, -1 long-way
A  = dm * sqrt(R1 * R2 * (1 + cos_dnu));

if abs(A) < 1e-8
    error('lambertSolver: transfer angle is 0 or 180 deg — degenerate case.');
end

% ---- scan z to find bracket ----
% For the 0-revolution case z in (-inf, (2*pi)^2).
% T(z) is monotonically decreasing, so there is exactly one crossing.
z_lo = -100.0;
z_hi = (2*pi)^2 - 1e-6;

N_scan = 200;
z_scan = linspace(z_lo, z_hi, N_scan);
t_scan = nan(1, N_scan);
for k = 1:N_scan
    t_scan(k) = tof_at_z(z_scan(k), R1, R2, A, mu);
end

% Find first valid bracket where T crosses desired tof
z_a = NaN;  z_b = NaN;
for k = 1:N_scan-1
    ta = t_scan(k);  tb = t_scan(k+1);
    if isfinite(ta) && isfinite(tb)
        if (ta - tof) * (tb - tof) <= 0
            z_a = z_scan(k);
            z_b = z_scan(k+1);
            break;
        end
    end
end

% If not found in positive-z range, widen search into deep hyperbolic region
if isnan(z_a)
    z_a = -1e4;
    z_b = z_scan(find(isfinite(t_scan), 1, 'first'));
    if isnan(z_b) || tof_at_z(z_a, R1, R2, A, mu) < tof
        error('lambertSolver: no solution found for TOF = %.1f s. Check inputs.', tof);
    end
end

% ---- bisection (60 iterations -> ~2^-60 relative error) ----
for k = 1:60
    z_mid = (z_a + z_b) / 2;
    t_mid = tof_at_z(z_mid, R1, R2, A, mu);
    if abs(t_mid - tof) / tof < 1e-10
        break;
    end
    if (tof_at_z(z_a, R1, R2, A, mu) - tof) * (t_mid - tof) <= 0
        z_b = z_mid;
    else
        z_a = z_mid;
    end
end
z_sol = (z_a + z_b) / 2;

% ---- Lagrange coefficients -> velocities ----
y_sol = y_fn(z_sol, R1, R2, A);
f     =  1 - y_sol / R1;
g     =  A * sqrt(y_sol / mu);
gdot  =  1 - y_sol / R2;

if abs(g) < 1e-12
    error('lambertSolver: degenerate g (near-singular geometry).');
end

v1 = (r2_vec - f * r1_vec) / g;
v2 = (gdot * r2_vec - r1_vec) / g;
end

% =========================================================
%  Local helpers
% =========================================================

function t = tof_at_z(z, R1, R2, A, mu)
    yv = y_fn(z, R1, R2, A);
    if yv <= 0
        t = inf;
        return;
    end
    Cv = stumpC(z);
    Sv = stumpS(z);
    xv = sqrt(yv / Cv);
    t  = (xv^3 * Sv + A * sqrt(yv)) / sqrt(mu);
    if ~isfinite(t) || t < 0
        t = inf;
    end
end

function y = y_fn(z, R1, R2, A)
    y = R1 + R2 + A * (z * stumpS(z) - 1) / sqrt(stumpC(z));
end

function C = stumpC(z)
    if z > 1e-6
        C = (1 - cos(sqrt(z))) / z;
    elseif z < -1e-6
        C = (cosh(sqrt(-z)) - 1) / (-z);
    else
        C = 0.5 - z/24 + z^2/720;
    end
end

function S = stumpS(z)
    if z > 1e-6
        sq = sqrt(z);
        S  = (sq - sin(sq)) / sq^3;
    elseif z < -1e-6
        sq = sqrt(-z);
        S  = (sinh(sq) - sq) / sq^3;
    else
        S = 1/6 - z/120 + z^2/5040;
    end
end
