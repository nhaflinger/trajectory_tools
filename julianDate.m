function jd = julianDate(year, month, day, hour, minute, second)
%JULIANDATE Convert calendar date/time to Julian Date (UTC)
%   jd = julianDate(year, month, day, hour, minute, second)
%   All inputs are scalars. Month is 1-12.

if nargin < 4 || isempty(hour), hour = 0; end
if nargin < 5 || isempty(minute), minute = 0; end
if nargin < 6 || isempty(second), second = 0; end

if month <= 2
    year = year - 1;
    month = month + 12;
end

A = floor(year/100);
B = 2 - A + floor(A/4);

jd = floor(365.25*(year+4716)) + floor(30.6001*(month+1)) + day + B - 1524.5;
jd = jd + (hour + minute/60 + second/3600) / 24;
end
