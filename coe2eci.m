function [r_vec, v_vec] = coe2eci(a, e, i, RAAN, omega, nu)
%COE2ECI  Classical orbital elements to ECI state vector (Earth).
%
%   [r_vec, v_vec] = coe2eci(a, e, i, RAAN, omega, nu)
%
%   Inputs (a in km, all angles in degrees):
%     a     - semi-major axis (km)
%     e     - eccentricity
%     i     - inclination (deg)
%     RAAN  - right ascension of ascending node (deg)
%     omega - argument of perigee (deg)
%     nu    - true anomaly (deg)
%
%   Outputs:
%     r_vec - ECI position (km, 3x1)
%     v_vec - ECI velocity (km/s, 3x1)
%
%   Uses Earth gravitational parameter mu = 398600.4418 km^3/s^2.

mu_E = 398600.4418;   % km^3/s^2

i_r  = deg2rad(i);
O_r  = deg2rad(RAAN);
om_r = deg2rad(omega);
nu_r = deg2rad(nu);

p = a * (1 - e^2);
r = p / (1 + e * cos(nu_r));

r_pqw = [r*cos(nu_r); r*sin(nu_r); 0];
v_pqw = sqrt(mu_E/p) * [-sin(nu_r); e + cos(nu_r); 0];

% Rotation: perifocal (PQW) -> ECI (same DCM as orbitalState.m)
cO = cos(O_r);  sO = sin(O_r);
co = cos(om_r); so = sin(om_r);
ci = cos(i_r);  si = sin(i_r);

R = [ cO*co - sO*so*ci,  -cO*so - sO*co*ci,  sO*si;
      sO*co + cO*so*ci,  -sO*so + cO*co*ci,  -cO*si;
      si*so,              si*co,               ci   ];

r_vec = R * r_pqw;
v_vec = R * v_pqw;
end
