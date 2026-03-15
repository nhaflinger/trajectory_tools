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
SSBodies.Earth = struct('name', 'Earth', ...
                        'mu', 398600.4418, ... % km^3/s^2
                        'radius', 6378.1363, ... % km (equatorial)
                        'a', 1.00000011 * SSBodies.Constants.AU, ... % semi-major axis of Earth's orbit (km)
                        'e', 0.01671022, ...
                        'inclination', 0.00005, ... % deg w.r.t. ecliptic
                        'obliquity', 23.45); % deg axial tilt to ecliptic

SSBodies.Moon = struct('name', 'Moon', ...
                        'mu', 4902.800066, ... % km^3/s^2
                        'radius', 1737.4, ... % km mean
                        'a', 384399, ... % km mean Earth-Moon distance
                        'e', 0.0549, ...
                        'inclination', 5.145, ... % degrees w.r.t. ecliptic
                        'obliquity', 6.68); % deg axial tilt to lunar orbital plane

%% Planets (simplified circular, coplanar orbits)
SSBodies.Mars = struct('name', 'Mars', ...
                        'mu', 4.282837e4, ... % km^3/s^2
                        'radius', 3389.5, ... % km mean
                        'a', 1.523679 * SSBodies.Constants.AU, ...
                        'e', 0.0934, ...
                        'inclination', 1.850, ...
                        'obliquity', 25.19); % deg axial tilt to ecliptic

SSBodies.Venus = struct('name', 'Venus', ...
                         'mu', 3.24859e5, ... % km^3/s^2
                         'radius', 6051.8, ... % km
                         'a', 0.723332 * SSBodies.Constants.AU, ...
                         'e', 0.0068, ...
                         'inclination', 3.394, ...
                         'obliquity', 177.4); % deg (retrograde rotation)

SSBodies.Jupiter = struct('name', 'Jupiter', ...
                          'mu', 1.26686534e8, ... % km^3/s^2
                          'radius', 69911, ... % km (equatorial)
                          'a', 5.204267 * SSBodies.Constants.AU, ...
                          'e', 0.0489, ...
                          'inclination', 1.305, ...
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
