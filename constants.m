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
SSBodies.Mercury = struct('name', 'Mercury', ...
                           'mu', 22031.86, ... % km^3/s^2
                           'radius', 2439.7, ... % km mean
                           'a', 0.38709893 * SSBodies.Constants.AU, ...
                           'e', 0.20563069, ...
                           'inclination', 7.00487, ... % deg w.r.t. ecliptic
                           'Omega',       48.33167, ... % deg, J2000
                           'omega_peri',  29.12478, ... % deg, J2000
                           'M0',         174.79439, ... % deg, J2000
                           'obliquity', 0.034);

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

%% Outer planets — J2000 Keplerian elements (Standish 1992, valid 1800-2050)
SSBodies.Saturn = struct('name', 'Saturn', ...
                         'mu', 3.793119e7, ... % km^3/s^2
                         'radius', 58232, ... % km (equatorial)
                         'a', 9.53707032 * SSBodies.Constants.AU, ...
                         'e', 0.05415060, ...
                         'inclination', 2.49117, ... % deg
                         'Omega',      113.71504, ... % deg, J2000
                         'omega_peri', 338.71690, ... % deg, J2000
                         'M0',         317.51238, ... % deg, J2000
                         'obliquity', 26.73);

SSBodies.Uranus = struct('name', 'Uranus', ...
                         'mu', 5.793832e6, ... % km^3/s^2
                         'radius', 25362, ... % km (equatorial)
                         'a', 19.19126393 * SSBodies.Constants.AU, ...
                         'e', 0.04716771, ...
                         'inclination', 0.76986, ... % deg
                         'Omega',       74.22988, ... % deg, J2000
                         'omega_peri',  96.73436, ... % deg, J2000
                         'M0',         142.26794, ... % deg, J2000
                         'obliquity', 97.77);

SSBodies.Neptune = struct('name', 'Neptune', ...
                          'mu', 6.836530e6, ... % km^3/s^2
                          'radius', 24622, ... % km (equatorial)
                          'a', 30.06896348 * SSBodies.Constants.AU, ...
                          'e', 0.00858587, ...
                          'inclination', 1.76917, ... % deg
                          'Omega',      131.72169, ... % deg, J2000
                          'omega_peri', 273.24966, ... % deg, J2000
                          'M0',         259.90868, ... % deg, J2000
                          'obliquity', 28.32);

SSBodies.Pluto = struct('name', 'Pluto', ...
                        'mu', 869.6, ... % km^3/s^2
                        'radius', 1188.3, ... % km mean
                        'a', 39.48168677 * SSBodies.Constants.AU, ...
                        'e', 0.24880766, ...
                        'inclination', 17.14175, ... % deg — highly inclined
                        'Omega',      110.30347, ... % deg, J2000
                        'omega_peri', 113.76329, ... % deg, J2000
                        'M0',          14.86205, ... % deg, J2000
                        'obliquity', 122.5);

%% Major asteroid belt bodies — J2000 elements (JPL Small Body Database)
SSBodies.Ceres = struct('name', 'Ceres', ...
                        'mu', 62.68, ... % km^3/s^2
                        'radius', 469.7, ... % km mean
                        'a', 2.76717 * SSBodies.Constants.AU, ...
                        'e', 0.07570, ...
                        'inclination', 10.5935, ... % deg
                        'Omega',       80.3268, ... % deg, J2000
                        'omega_peri',  73.5976, ... % deg, J2000
                        'M0',          80.3929, ... % deg, J2000
                        'obliquity', 4.0);

SSBodies.Vesta = struct('name', 'Vesta', ...
                        'mu', 17.288, ... % km^3/s^2
                        'radius', 262.7, ... % km mean
                        'a', 2.36179 * SSBodies.Constants.AU, ...
                        'e', 0.08874, ...
                        'inclination', 7.1401, ... % deg
                        'Omega',      103.8514, ... % deg, J2000
                        'omega_peri', 151.1985, ... % deg, J2000
                        'M0',          20.8632, ... % deg, J2000
                        'obliquity', 29.1);

SSBodies.Pallas = struct('name', 'Pallas', ...
                         'mu', 14.0, ... % km^3/s^2
                         'radius', 272.5, ... % km mean
                         'a', 2.77231 * SSBodies.Constants.AU, ...
                         'e', 0.23048, ...
                         'inclination', 34.8355, ... % deg — highly inclined
                         'Omega',      173.0962, ... % deg, J2000
                         'omega_peri', 310.0521, ... % deg, J2000
                         'M0',          78.2272, ... % deg, J2000
                         'obliquity', 84.0);

SSBodies.Hygiea = struct('name', 'Hygiea', ...
                         'mu', 5.628, ... % km^3/s^2
                         'radius', 217.0, ... % km mean
                         'a', 3.13927 * SSBodies.Constants.AU, ...
                         'e', 0.11754, ...
                         'inclination', 3.8317, ... % deg
                         'Omega',      283.4956, ... % deg, J2000
                         'omega_peri', 312.3112, ... % deg, J2000
                         'M0',         178.2601, ... % deg, J2000
                         'obliquity', 0.0);

%% Trans-Neptunian dwarf planets (approximate J2000 elements, JPL Horizons)
SSBodies.Eris = struct('name', 'Eris', ...
                       'mu', 1108.0, ... % km^3/s^2
                       'radius', 1163.0, ... % km mean
                       'a', 67.781 * SSBodies.Constants.AU, ...
                       'e', 0.44068, ...
                       'inclination', 44.0445, ... % deg — highly inclined
                       'Omega',       35.9531, ... % deg, J2000
                       'omega_peri', 151.4305, ... % deg, J2000
                       'M0',         205.9823, ... % deg, J2000
                       'obliquity', 78.0);

SSBodies.Makemake = struct('name', 'Makemake', ...
                           'mu', 300.0, ... % km^3/s^2 (uncertain)
                           'radius', 715.0, ... % km mean (approximate)
                           'a', 45.430 * SSBodies.Constants.AU, ...
                           'e', 0.16254, ...
                           'inclination', 28.9834, ... % deg
                           'Omega',       79.3551, ... % deg, J2000
                           'omega_peri', 295.6536, ... % deg, J2000
                           'M0',         168.2614, ... % deg, J2000
                           'obliquity', 0.0);

SSBodies.Haumea = struct('name', 'Haumea', ...
                         'mu', 330.0, ... % km^3/s^2 (uncertain)
                         'radius', 780.0, ... % km mean (triaxial; use mean)
                         'a', 43.116 * SSBodies.Constants.AU, ...
                         'e', 0.19126, ...
                         'inclination', 28.2137, ... % deg
                         'Omega',      122.1670, ... % deg, J2000
                         'omega_peri', 239.1781, ... % deg, J2000
                         'M0',         218.2656, ... % deg, J2000
                         'obliquity', 0.0);

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
