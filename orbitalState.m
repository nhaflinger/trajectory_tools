function [r_vec, v_vec] = orbitalState(body, jd, epoch0)
%ORBITALSTATE  Heliocentric ecliptic state vector from Keplerian elements.
%
%   [r_vec, v_vec] = orbitalState(body, jd)
%   [r_vec, v_vec] = orbitalState(body, jd, epoch0)
%
%   Propagates body's J2000 Keplerian elements to Julian Date jd and
%   returns position (km) and velocity (km/s) in the heliocentric ecliptic
%   J2000 frame.  No Aerospace Toolbox required.
%
%   body must contain:
%     a          - semi-major axis (km)
%     e          - eccentricity
%     inclination- orbital inclination (deg)
%     Omega      - longitude of ascending node (deg), J2000
%     omega_peri - argument of perihelion (deg), J2000
%     M0         - mean anomaly at epoch0 (deg)
%
%   If Omega / omega_peri / M0 are absent the function falls back to a
%   circular, coplanar approximation (same as the old orbitalStateCircular).

if nargin < 3 || isempty(epoch0)
    epoch0 = 2451545.0;   % J2000.0
end

muSun = constants().Sun.mu;   % km^3/s^2

% ---- hyperbolic orbit (e >= 1) ----
% Requires body.t_peri_jd  (Julian Date of perihelion passage)
%          body.e, body.inclination, body.Omega, body.omega_peri
% body.a may be negative (convention) or positive; abs value is used.
if isfield(body,'e') && body.e >= 1.0
    abs_a = abs(body.a);                         % |semi-major axis|
    n_h   = sqrt(muSun / abs_a^3);               % hyperbolic mean motion
    dt    = (jd - body.t_peri_jd) * 86400;       % seconds from perihelion (signed)
    M_h   = n_h * dt;                            % hyperbolic mean anomaly (signed)
    H     = hyperbolicSolve(M_h, body.e);        % hyperbolic anomaly
    nu    = 2 * atan(sqrt((body.e+1)/(body.e-1)) * tanh(H/2));
    p     = abs_a * (body.e^2 - 1);
    r     = p / (1 + body.e * cos(nu));
    r_pqw = r * [cos(nu); sin(nu); 0];
    v_pqw = sqrt(muSun/p) * [-sin(nu); body.e + cos(nu); 0];
    % Rotate perifocal -> heliocentric ecliptic (same DCM as elliptic)
    Om  = deg2rad(body.Omega);
    om  = deg2rad(body.omega_peri);
    inc = deg2rad(body.inclination);
    cOm = cos(Om);  sOm = sin(Om);
    com = cos(om);  som = sin(om);
    ci  = cos(inc); si  = sin(inc);
    R = [ cOm*com - sOm*som*ci,  -cOm*som - sOm*com*ci,  sOm*si;
          sOm*com + cOm*som*ci,  -sOm*som + cOm*com*ci,  -cOm*si;
          si*som,                  si*com,                  ci    ];
    r_vec = R * r_pqw;
    v_vec = R * v_pqw;
    return;
end

% ---- fallback: circular coplanar ----
if ~isfield(body,'Omega') || ~isfield(body,'omega_peri') || ~isfield(body,'M0')
    n     = sqrt(muSun / body.a^3);
    dt    = (jd - epoch0) * 86400;
    theta = mod(n*dt, 2*pi);
    r_vec = [body.a*cos(theta); body.a*sin(theta); 0];
    vMag  = sqrt(muSun / body.a);
    v_vec = vMag * [-sin(theta); cos(theta); 0];
    return;
end

% ---- propagate Keplerian elements ----
a  = body.a;
e  = body.e;
n  = sqrt(muSun / a^3);              % mean motion (rad/s)
dt = (jd - epoch0) * 86400;          % seconds since epoch

M  = deg2rad(body.M0) + n*dt;
M  = mod(M, 2*pi);

% Kepler's equation M = E - e*sin(E) — Newton-Raphson
E = keplerSolve(M, e);

% True anomaly
sinE = sin(E);  cosE = cos(E);
nu   = atan2(sqrt(1 - e^2) * sinE, cosE - e);

% Radius and specific angular momentum
r = a * (1 - e*cosE);
p = a * (1 - e^2);
h = sqrt(muSun * p);

% State in perifocal (PQW) frame: x along periapsis, y 90 deg ahead
r_pqw = r  * [cos(nu); sin(nu); 0];
v_pqw = (muSun/h) * [-sin(nu); e + cos(nu); 0];

% Rotation: perifocal -> heliocentric ecliptic.
% Uses the direction cosine matrix from Bate, Mueller & White (1971)
% eqs. 4.6-5, which correctly maps PQW -> ECI/ecliptic:
%
%   R(i,j) = cosine of angle between i-th ecliptic base vector
%             and j-th perifocal base vector.
%
Om  = deg2rad(body.Omega);
om  = deg2rad(body.omega_peri);
inc = deg2rad(body.inclination);

cOm = cos(Om);  sOm = sin(Om);
com = cos(om);  som = sin(om);
ci  = cos(inc); si  = sin(inc);

R = [ cOm*com - sOm*som*ci,  -cOm*som - sOm*com*ci,  sOm*si;
      sOm*com + cOm*som*ci,  -sOm*som + cOm*com*ci,  -cOm*si;
      si*som,                  si*com,                  ci    ];

r_vec = R * r_pqw;
v_vec = R * v_pqw;
end

% ---- local helper ----

function E = keplerSolve(M, e)
%KEPLERSOLVE  Newton-Raphson solution of Kepler's equation M = E - e*sin(E).
    E = M;
    for k = 1:50
        dE = (M - E + e*sin(E)) / (1 - e*cos(E));
        E  = E + dE;
        if abs(dE) < 1e-13, break; end
    end
end

function H = hyperbolicSolve(M_h, e)
%HYPERBOLICSOLVE  Newton-Raphson solution of hyperbolic Kepler's equation
%   M_h = e*sinh(H) - H.  Works for any signed M_h.
    H = sign(M_h) * log(2*abs(M_h)/e + 1.8);   % Battin initial guess
    for k = 1:50
        dH = (M_h - e*sinh(H) + H) / (e*cosh(H) - 1);
        H  = H + dH;
        if abs(dH) < 1e-13, break; end
    end
end
