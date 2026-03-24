function coe = eci2coe(r_vec, v_vec)
%ECI2COE  ECI state vector to classical orbital elements (Earth orbit).
%
%   coe = eci2coe(r_vec, v_vec)
%
%   Inputs:
%     r_vec - ECI position (km, 3x1 or 1x3)
%     v_vec - ECI velocity (km/s, 3x1 or 1x3)
%
%   Output struct fields (angles in degrees, distances in km):
%     a     - semi-major axis (km)
%     e     - eccentricity
%     i     - inclination (deg)
%     RAAN  - right ascension of ascending node (deg)
%     omega - argument of perigee (deg)
%     nu    - true anomaly (deg)
%     M     - mean anomaly (deg)
%     p     - semi-latus rectum (km)
%     h_vec - specific angular momentum vector (km^2/s, 3x1)

mu_E = 398600.4418;   % km^3/s^2

r_vec = r_vec(:);
v_vec = v_vec(:);
r = norm(r_vec);
v = norm(v_vec);

% Specific angular momentum
h_vec = cross(r_vec, v_vec);
h     = norm(h_vec);

% Eccentricity vector (points toward periapsis)
e_vec = cross(v_vec, h_vec)/mu_E - r_vec/r;
e     = norm(e_vec);

% Node vector (K x h, points toward ascending node)
K     = [0; 0; 1];
N_vec = cross(K, h_vec);
N     = norm(N_vec);

% Inclination
i = acosd(max(-1, min(1, h_vec(3)/h)));

% RAAN
if N > 1e-10
    RAAN = acosd(max(-1, min(1, N_vec(1)/N)));
    if N_vec(2) < 0, RAAN = 360 - RAAN; end
else
    RAAN = 0;   % equatorial orbit: undefined, set to 0
end

% Argument of perigee
if e > 1e-10
    if N > 1e-10
        omega = acosd(max(-1, min(1, dot(N_vec, e_vec)/(N*e))));
        if e_vec(3) < 0, omega = 360 - omega; end
    else
        % equatorial elliptic
        omega = atan2d(e_vec(2), e_vec(1));
        if omega < 0, omega = omega + 360; end
    end
else
    omega = 0;  % circular: undefined, set to 0
end

% True anomaly
if e > 1e-10
    nu = acosd(max(-1, min(1, dot(e_vec, r_vec)/(e*r))));
    if dot(r_vec, v_vec) < 0, nu = 360 - nu; end
else
    % circular: argument of latitude from node
    if N > 1e-10
        nu = acosd(max(-1, min(1, dot(N_vec, r_vec)/(N*r))));
        if r_vec(3) < 0, nu = 360 - nu; end
    else
        nu = atan2d(r_vec(2), r_vec(1));
        if nu < 0, nu = nu + 360; end
    end
end

% Semi-major axis (vis-viva equation)
a = 1 / (2/r - v^2/mu_E);

% Semi-latus rectum
p = h^2 / mu_E;

% Mean anomaly (through eccentric anomaly)
nu_r = deg2rad(nu);
cosE = (e + cos(nu_r)) / (1 + e*cos(nu_r));
sinE = sqrt(max(0, 1 - e^2)) * sin(nu_r) / (1 + e*cos(nu_r));
E    = atan2(sinE, cosE);
M    = rad2deg(mod(E - e*sin(E), 2*pi));

coe = struct('a', a, 'e', e, 'i', i, 'RAAN', RAAN, 'omega', omega, ...
             'nu', nu, 'M', M, 'p', p, 'h_vec', h_vec);
end
