function [r, v] = orbitalStateCircular(body, epoch, epoch0)
%ORBITALSTATECIRCULAR  Heliocentric state for a body around the Sun.
%
%   Delegates to orbitalState when the body struct contains full J2000
%   Keplerian elements (Omega, omega_peri, M0), giving an eccentric 3-D
%   solution.  Falls back to the original circular-coplanar approximation
%   for bodies that only carry semi-major axis data.

if nargin < 3 || isempty(epoch0)
    epoch0 = 2451545.0;  % J2000
end

[r, v] = orbitalState(body, epoch, epoch0);
end
