function orb = earthOrbit(type, varargin)
%EARTHORBIT  Create an Earth orbit definition struct.
%
%   orb = earthOrbit('coe', a_km, e, i_deg, RAAN_deg, omega_deg, M0_deg)
%   orb = earthOrbit('coe', a_km, e, i_deg, RAAN_deg, omega_deg, M0_deg, epoch_jd)
%   orb = earthOrbit('eci', r_vec, v_vec)
%   orb = earthOrbit('eci', r_vec, v_vec, epoch_jd)
%   orb = earthOrbit('circular', alt_km, i_deg)
%   orb = earthOrbit('circular', alt_km, i_deg, RAAN_deg)
%   orb = earthOrbit('sso', alt_km)
%   orb = earthOrbit('sso', alt_km, e, omega_deg)
%   orb = earthOrbit('sso', alt_km, e, omega_deg, ltan_hrs, epoch_jd)
%   orb = earthOrbit('geo')
%   orb = earthOrbit('geo', lon_deg)
%   orb = earthOrbit('molniya', RAAN_deg)
%   orb = earthOrbit('molniya', RAAN_deg, omega_deg)
%   orb = earthOrbit('tle', line1, line2)
%
%   SSO LTAN note: ltan_hrs is the Local Time of Ascending Node in hours (0-24).
%   When provided, the RAAN is set so that the ascending node aligns with the
%   specified local solar time at the given epoch_jd. epoch_jd defaults to J2000
%   if not specified, which means the RAAN will be correct for launches near
%   Jan 1, 2000; for other launch dates, pass the actual epoch JD.
%
%   Returns a struct with fields:
%     a              - semi-major axis (km)
%     e              - eccentricity
%     i              - inclination (deg)
%     RAAN           - right ascension of ascending node (deg)
%     omega          - argument of perigee (deg)
%     M0             - mean anomaly at epoch (deg)
%     epoch_jd       - epoch as Julian Date
%     period         - orbital period (s)
%     alt_peri       - periapsis altitude (km)
%     alt_apo        - apoapsis altitude (km)
%     r_vec          - ECI position at epoch (km, 3x1)
%     v_vec          - ECI velocity at epoch (km/s, 3x1)
%     type           - orbit type label
%     ltan_hrs       - LTAN in hours (NaN if not specified / non-SSO)
%     sun_sync_RAAN  - RAAN set from LTAN calculation (NaN if not applicable)

mu_E  = 398600.4418;    % km^3/s^2
R_E   = 6378.1363;      % km
J2    = 1.08262668e-3;

% Default LTAN fields (set to NaN; overridden by SSO+LTAN case below)
ltan_hrs      = NaN;
sun_sync_RAAN = NaN;

type = lower(strtrim(type));

switch type

    %% ── Classical orbital elements ─────────────────────────────────────────
    case 'coe'
        if numel(varargin) < 6
            error('earthOrbit: ''coe'' requires a, e, i, RAAN, omega, M0');
        end
        a     = varargin{1};
        e     = varargin{2};
        i     = varargin{3};
        RAAN  = varargin{4};
        omega = varargin{5};
        M0    = varargin{6};
        epoch_jd = getEpoch(varargin, 7);

    %% ── ECI state vector ───────────────────────────────────────────────────
    case 'eci'
        if numel(varargin) < 2
            error('earthOrbit: ''eci'' requires r_vec, v_vec');
        end
        coe_s = eci2coe(varargin{1}(:), varargin{2}(:));
        a     = coe_s.a;
        e     = coe_s.e;
        i     = coe_s.i;
        RAAN  = coe_s.RAAN;
        omega = coe_s.omega;
        M0    = coe_s.M;
        epoch_jd = getEpoch(varargin, 3);

    %% ── Circular orbit at given altitude and inclination ───────────────────
    case 'circular'
        if numel(varargin) < 2
            error('earthOrbit: ''circular'' requires alt_km, i_deg');
        end
        a     = R_E + varargin{1};
        e     = 0;
        i     = varargin{2};
        RAAN  = 0;   if numel(varargin) >= 3, RAAN = varargin{3}; end
        omega = 0;
        M0    = 0;
        epoch_jd = 2451545.0;

    %% ── Sun-synchronous orbit ──────────────────────────────────────────────
    case 'sso'
        if numel(varargin) < 1
            error('earthOrbit: ''sso'' requires alt_km');
        end
        alt   = varargin{1};
        e     = 0;  if numel(varargin) >= 2, e     = varargin{2}; end
        omega = 0;  if numel(varargin) >= 3, omega = varargin{3}; end
        a     = R_E + alt;
        RAAN  = 0;
        M0    = 0;
        % Required nodal regression = Earth's mean orbital rate around Sun
        om_sun = 2*pi / (365.2422 * 86400);     % rad/s
        p_orb  = a * (1 - e^2);
        n      = sqrt(mu_E / a^3);
        cos_i  = -2/3 * (p_orb/R_E)^2 * om_sun / (n * J2);
        if abs(cos_i) > 1
            error('earthOrbit: SSO not achievable at alt = %.0f km (cos_i = %.3f out of range)', alt, cos_i);
        end
        i = acosd(cos_i);
        epoch_jd = 2451545.0;

        % LTAN (Local Time of Ascending Node) support
        if numel(varargin) >= 4
            ltan_hrs = varargin{4};
            if numel(varargin) >= 5
                epoch_jd = varargin{5};
            end
            % Sun RA at epoch
            [sun_hat, ~] = sunPosition(epoch_jd);
            sun_RA = atan2d(sun_hat(2), sun_hat(1));   % deg, -180..180
            % Required RAAN: LTAN=12h means ascending node at noon (Sun's RA)
            RAAN = mod(sun_RA + (ltan_hrs - 12) * 15, 360);
            sun_sync_RAAN = RAAN;
        end

    %% ── Geostationary orbit ────────────────────────────────────────────────
    case 'geo'
        om_E  = 7.2921150e-5;              % rad/s Earth rotation
        a     = (mu_E / om_E^2)^(1/3);    % 42164.17 km
        e     = 0;
        i     = 0;
        omega = 0;
        M0    = 0;
        RAAN  = 0;
        if numel(varargin) >= 1
            % Place satellite over given geodetic longitude at epoch
            % For i=0, e=0: satellite longitude = RAAN + M0 - GMST_epoch.
            % Simplification: set RAAN = lon_deg (valid near J2000 GMST = 280.46 deg).
            RAAN = varargin{1};
        end
        epoch_jd = 2451545.0;

    %% ── Molniya orbit ──────────────────────────────────────────────────────
    case 'molniya'
        % Standard Molniya: i=63.435°, e=0.74, T=12 h, omega=270° (apogee over NH)
        T_mol = 12 * 3600;                          % 12-hour period (s)
        a     = (mu_E * (T_mol / (2*pi))^2)^(1/3);  % ~26560 km
        e     = 0.74;
        i     = 63.435;
        omega = 270;
        M0    = 0;
        RAAN  = 0;  if numel(varargin) >= 1, RAAN  = varargin{1}; end
        omega = 270; if numel(varargin) >= 2, omega = varargin{2}; end
        epoch_jd = 2451545.0;

    %% ── Two-Line Elements ──────────────────────────────────────────────────
    case 'tle'
        if numel(varargin) < 2
            error('earthOrbit: ''tle'' requires line1, line2');
        end
        [a, e, i, RAAN, omega, M0, epoch_jd] = parseTLE(varargin{1}, varargin{2});

    otherwise
        error('earthOrbit: unknown type ''%s''. Valid: coe, eci, circular, sso, geo, molniya, tle', type);
end

%% ── Derived quantities ─────────────────────────────────────────────────────
period   = 2*pi * sqrt(a^3 / mu_E);   % seconds
rp       = a * (1 - e);
ra       = a * (1 + e);
alt_peri = rp - R_E;
alt_apo  = ra - R_E;

% ECI state at epoch: convert M0 -> E -> nu -> state
E0  = keplerSolve(deg2rad(M0), e);
nu0 = 2 * atan2(sqrt(1+e) * sin(E0/2), sqrt(1-e) * cos(E0/2));
[r_vec, v_vec] = coe2eci(a, e, i, RAAN, omega, rad2deg(nu0));

orb = struct( ...
    'a',             a,             ...
    'e',             e,             ...
    'i',             i,             ...
    'RAAN',          RAAN,          ...
    'omega',         omega,         ...
    'M0',            M0,            ...
    'epoch_jd',      epoch_jd,      ...
    'period',        period,        ...
    'alt_peri',      alt_peri,      ...
    'alt_apo',       alt_apo,       ...
    'r_vec',         r_vec,         ...
    'v_vec',         v_vec,         ...
    'type',          type,          ...
    'ltan_hrs',      ltan_hrs,      ...
    'sun_sync_RAAN', sun_sync_RAAN);
end

%% ── Local helpers ──────────────────────────────────────────────────────────

function epoch_jd = getEpoch(args, idx)
    if numel(args) >= idx
        epoch_jd = args{idx};
    else
        epoch_jd = 2451545.0;   % J2000.0
    end
end

function E = keplerSolve(M, e)
    E = M;
    for k = 1:50
        dE = (M - E + e*sin(E)) / (1 - e*cos(E));
        E  = E + dE;
        if abs(dE) < 1e-13, break; end
    end
end

function [a, e, i, RAAN, omega, M0, epoch_jd] = parseTLE(line1, line2)
%PARSETLE  Extract mean elements from a two-line element set (SGP4 mean elements).
%   Column indices follow the standard TLE format as documented by CelesTrak.
    mu_E = 398600.4418;

    % Epoch from line 1 (cols 19-20: 2-digit year, cols 21-32: fractional day)
    yr2      = str2double(line1(19:20));
    yr       = 2000 + yr2;
    if yr2 >= 57, yr = 1900 + yr2; end   % pre-1957 orbits won't exist
    doy      = str2double(line1(21:32));  % fractional day of year
    epoch_jd = julianDate(yr, 1, 1) + (doy - 1);

    % Elements from line 2
    i     = str2double(line2(9:16));
    RAAN  = str2double(line2(18:25));
    e     = str2double(['0.' strtrim(line2(27:33))]);
    omega = str2double(line2(35:42));
    M0    = str2double(line2(44:51));
    n_rev = str2double(line2(53:63));    % rev/day
    n     = n_rev * 2*pi / 86400;       % rad/s
    a     = (mu_E / n^2)^(1/3);         % km
end
