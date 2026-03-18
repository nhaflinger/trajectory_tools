function SSBodies = constants()
%CONSTANTS Solar system body constants for patched conic approximations
%   SSBodies is a struct containing physical and orbital constants for
%   planets, the Moon, and the Sun. Units are SI unless otherwise noted.
%
%   Usage:
%       bodies = constants();
%       muEarth = bodies.Earth.mu;
%       r_SOI = bodies.Moon.soi;

%% Universal constants
SSBodies.Constants = struct('G', 6.67430e-20, ... % km^3 / kg / s^2
                            'AU', 149597870.7, ... % km
                            'day', 86400); % s

%% Sun
SSBodies.Sun = struct('name', 'Sun', ...
                      'mu', 1.32712440018e11, ... % km^3/s^2
                      'radius', 695700); % km (photosphere)

%% Earth-Moon system
% J2000 Keplerian elements from Standish (1992), valid 1800-2050.
% Omega, omega_peri, M0 derived from: L0 (mean longitude) and
% omega_bar (longitude of perihelion) tabulated by JPL:
%   omega_peri = omega_bar - Omega
%   M0         = L0 - omega_bar   (wrapped to [0,360))
SSBodies.Earth = struct('name', 'Earth', ...
                        'mu', 398600.4418, ... % km^3/s^2
                        'radius', 6378.1363, ... % km (equatorial)
                        'a', 1.00000011 * SSBodies.Constants.AU, ... % semi-major axis (km)
                        'e', 0.01671022, ...
                        'inclination', 0.00005, ... % deg w.r.t. ecliptic
                        'Omega',       -11.26064, ... % deg, longitude of ascending node, J2000
                        'omega_peri',  114.20783, ... % deg, argument of perihelion, J2000
                        'M0',          357.51716, ... % deg, mean anomaly at J2000.0
                        'obliquity', 23.45); % deg axial tilt to ecliptic

SSBodies.Moon = struct('name', 'Moon', ...
                        'mu', 4902.800066, ... % km^3/s^2
                        'radius', 1737.4, ... % km mean
                        'a', 384399, ... % km mean Earth-Moon distance
                        'e', 0.0549, ...
                        'inclination', 5.145, ... % degrees w.r.t. ecliptic
                        'obliquity', 6.68); % deg axial tilt to lunar orbital plane

%% Planets — J2000 Keplerian elements (Standish 1992, valid 1800-2050)
SSBodies.Mars = struct('name', 'Mars', ...
                        'mu', 4.282837e4, ... % km^3/s^2
                        'radius', 3389.5, ... % km mean
                        'a', 1.52366231 * SSBodies.Constants.AU, ...
                        'e', 0.09341233, ...
                        'inclination', 1.85061, ... % deg w.r.t. ecliptic
                        'Omega',       49.57854, ... % deg, J2000
                        'omega_peri', 286.46230, ... % deg, J2000
                        'M0',          19.41248, ... % deg, J2000
                        'obliquity', 25.19);

SSBodies.Venus = struct('name', 'Venus', ...
                         'mu', 3.24859e5, ... % km^3/s^2
                         'radius', 6051.8, ... % km
                         'a', 0.72332102 * SSBodies.Constants.AU, ...
                         'e', 0.00676399, ...
                         'inclination', 3.39777, ... % deg w.r.t. ecliptic
                         'Omega',       76.67069, ... % deg, J2000
                         'omega_peri',  54.86229, ... % deg, J2000
                         'M0',          50.44675, ... % deg, J2000
                         'obliquity', 177.4);

SSBodies.Jupiter = struct('name', 'Jupiter', ...
                          'mu', 1.26686534e8, ... % km^3/s^2
                          'radius', 69911, ... % km (equatorial)
                          'a', 5.20336301 * SSBodies.Constants.AU, ...
                          'e', 0.04839266, ...
                          'inclination', 1.30530, ... % deg w.r.t. ecliptic
                          'Omega',      100.55615, ... % deg, J2000
                          'omega_peri', 273.71880, ... % deg, J2000
                          'M0',          20.05984, ... % deg, J2000
                          'obliquity', 3.13);

%% Major asteroid belt bodies
SSBodies.Ceres = struct('name', 'Ceres', ...
                        'mu', 62.68, ... % km^3/s^2
                        'radius', 469.7, ... % km mean
                        'a', 2.7675 * SSBodies.Constants.AU, ...
                        'e', 0.0758, ...
                        'inclination', 10.593, ... % deg
                        'obliquity', 4.0); % deg

SSBodies.Vesta = struct('name', 'Vesta', ...
                        'mu', 17.288, ... % km^3/s^2
                        'radius', 262.7, ... % km mean
                        'a', 2.3615 * SSBodies.Constants.AU, ...
                        'e', 0.0887, ...
                        'inclination', 7.140, ...
                        'obliquity', 29.1);

SSBodies.Pallas = struct('name', 'Pallas', ...
                         'mu', 14.0, ... % km^3/s^2
                         'radius', 272.5, ... % km mean
                         'a', 2.7724 * SSBodies.Constants.AU, ...
                         'e', 0.2313, ...
                         'inclination', 34.841, ...
                         'obliquity', 84.0); % highly tilted

SSBodies.Hygiea = struct('name', 'Hygiea', ...
                         'mu', 5.628, ... % km^3/s^2
                         'radius', 217.0, ... % km mean
                         'a', 3.1417 * SSBodies.Constants.AU, ...
                         'e', 0.1125, ...
                         'inclination', 3.832, ...
                         'obliquity', 0.0); % unknown, assumed negligible

%% Sphere-of-influence helper
bodies = fieldnames(SSBodies);
for i = 1:numel(bodies)
    body = SSBodies.(bodies{i});
    if isfield(body, 'a') && strcmp(body.name, 'Sun') == false
        % Approximate SOI relative to Sun
        muPrimary = SSBodies.Sun.mu;
        body.soi = body.a * (body.mu / muPrimary)^(2/5);
        SSBodies.(bodies{i}) = body;
    end
end

end
