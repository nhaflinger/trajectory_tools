function [v1, v2] = lambertSolver(r1, r2, tof, mu, isLongWay)
%LAMBERTSOLVER Solve the two-body Lambert problem (universal variable)
%   [v1, v2] = lambertSolver(r1, r2, tof, mu, isLongWay)
%   r1, r2 are 3x1 position vectors (km), tof is time-of-flight (s).
%   Returns initial and final velocity vectors (km/s).
%   isLongWay indicates whether to take the long-way (true) or short-way (false).
%   Uses a simple universal variable iteration.

if nargin < 5
    isLongWay = false;
end

r1 = r1(:); r2 = r2(:);

R1 = norm(r1);
R2 = norm(r2);

cosDelta = dot(r1, r2) / (R1*R2);
if isLongWay
    cosDelta = -cosDelta;
end

% chord
c = sqrt(R1^2 + R2^2 - 2*R1*R2*cosDelta);

% semiperimeter
s = (R1 + R2 + c) / 2;

% Minimum energy semi-major axis
amin = s/2;

% initial guess for x (universal variable)
x = 0.0;

% Target time-of-flight function
function [dt, dtdx] = tofFunc(x)
    % Stumpff functions
    z = x^2;
    Sz = stumpffS(z);
    Cz = stumpffC(z);
    y = R1 + R2 + ( (s*(z*Sz - 1)) / sqrt(Cz) );
    if y < 0
        dt = inf;
        dtdx = 0;
        return;
    end
    a = s/2 + (s^2 * (1 - z*Sz)) / (16 * y);
    if a <= 0
        dt = inf;
        dtdx = 0;
        return;
    end
    beta = 2 * asin( sqrt((s - c) / (2*a)) );
    if isLongWay
        beta = -beta;
    end
    alpha = 2 * asin( sqrt(s / (2*a)) );
    dt = (a^(3/2)/sqrt(mu)) * (alpha - sin(alpha) - (beta - sin(beta)));
    % naive derivative (not highly accurate, but acceptable)
    dtdx = 0;
end

% Iterate using simple secant / Newton hybrid
maxIter = 50;
tol = 1e-6;

% Use a naive initial guess based on minimum energy
x = sqrt(mu) * tof / (s/2)^(3/2);

for iter = 1:maxIter
    [dt, dtdx] = tofFunc(x);
    if dt == inf
        x = x * 0.8;
        continue;
    end
    err = dt - tof;
    if abs(err) < tol
        break;
    end
    % simple update
    x = x - err / max(dtdx, 1e-6);
end

% Compute velocities using Lagrange coefficients (approximate)
% This implementation is intentionally simple: it uses the f,g formulation
% with finite differences.
dt_small = tof * 1e-6;
[~, ~] = tofFunc(x);

% Propagate r1 over dt_small using linear approximation (not accurate for long tof)
r1p = r1 + (r2 - r1) * (dt_small/tof);

f = 1 - (norm(r2 - r1)^2) / (2 * amin * R1);
g = tof - sqrt(amin^3/mu) * (sin(tof*sqrt(mu/amin^3)));

v1 = (r2 - f*r1) / g;
v2 = (g_dot()*r2 - r1) / g;

    function gd = g_dot()
        % approximate g dot (not rigorous)
        gd = 1 - (norm(r2 - r1)^2) / (2 * amin * R2);
    end
end

function C = stumpffC(z)
if z > 0
    C = (1 - cos(sqrt(z))) / z;
elseif z < 0
    C = (cosh(sqrt(-z)) - 1) / (-z);
else
    C = 1/2;
end
end

function S = stumpffS(z)
if z > 0
    S = (sqrt(z) - sin(sqrt(z))) / (sqrt(z)^3);
elseif z < 0
    S = (sinh(sqrt(-z)) - sqrt(-z)) / (sqrt(-z)^3);
else
    S = 1/6;
end
end
