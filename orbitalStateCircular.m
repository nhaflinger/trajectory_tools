function [r,v] = orbitalStateCircular(body, epoch, epoch0)
%ORBITALSTATECIRCULAR Circular orbital state for a body around the Sun
%   [r,v] = orbitalStateCircular(body, epoch, epoch0)
%   Returns a simple circular-coplanar heliocentric state vector (km, km/s).
%   epoch and epoch0 are in Julian days. epoch0 defines the reference
%   epoch at which true anomaly = 0.
%
%   body should be a struct with fields: a (semi-major axis, km), mu (km^3/s^2)
%   and optionally meanMotion (rad/s).

if nargin < 3 || isempty(epoch0)
    epoch0 = 2451545.0; % J2000
end

% Mean motion (rad/s)
if isfield(body, 'meanMotion')
    n = body.meanMotion;
else
    n = sqrt(constants().Sun.mu / body.a^3);
end

% Time since epoch0 in seconds
dt = (epoch - epoch0) * 86400;

% True anomaly (circular orbit)
theta = mod(n * dt, 2*pi);

% Position in ecliptic plane
r = [body.a * cos(theta); body.a * sin(theta); 0];

% Velocity magnitude
vMag = sqrt(constants().Sun.mu / body.a);

% Velocity perpendicular to position (prograde)
v = vMag * [-sin(theta); cos(theta); 0];
end
