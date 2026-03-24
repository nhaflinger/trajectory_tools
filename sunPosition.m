function [r_sun_hat, r_sun_km] = sunPosition(jd)
%SUNPOSITION  Compute the unit vector from Earth center to Sun in ECI (J2000).
%
%   [r_sun_hat, r_sun_km] = sunPosition(jd)
%
%   Low-precision solar almanac accurate to ~0.01 deg, sufficient for
%   beta angle and eclipse calculations.
%
%   Algorithm from Vallado "Fundamentals of Astrodynamics and Applications"
%   4th ed, Algorithm 29.
%
%   Input:
%     jd        - Julian Date (scalar)
%
%   Outputs:
%     r_sun_hat - unit vector from Earth to Sun in ECI J2000, 3x1
%     r_sun_km  - full position vector from Earth to Sun in ECI J2000 (km), 3x1

AU_km = 149597870.7;   % km per AU

T = (jd - 2451545.0) / 36525;   % Julian centuries from J2000

M_sun  = 357.5291092  + 35999.0502909 * T;   % mean anomaly (deg)
L_sun  = 280.4664567  + 36000.76982   * T;   % mean longitude (deg)

lambda = L_sun + 1.914666471 * sind(M_sun) + 0.019994643 * sind(2*M_sun);  % ecliptic longitude (deg)
eps    = 23.439291111 - 0.013004167 * T;     % obliquity of ecliptic (deg)
r_AU   = 1.000140612  - 0.016708617 * cosd(M_sun) - 0.000139589 * cosd(2*M_sun);  % Sun-Earth distance (AU)

r_sun = r_AU * AU_km * [cosd(lambda); ...
                         sind(lambda)*cosd(eps); ...
                         sind(lambda)*sind(eps)];

r_sun_hat = r_sun / norm(r_sun);
r_sun_km  = r_sun;
end
