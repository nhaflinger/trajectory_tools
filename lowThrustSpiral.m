function result = lowThrustSpiral(body, r0_km, r1_km, options)
%LOWTHRUSTSPIRAL  Edelbaum analytical spiral + tangential-thrust RK4 propagation.
%
%   result = lowThrustSpiral(body, r0_km, r1_km)
%   result = lowThrustSpiral(body, r0_km, r1_km, options)
%
%   Computes the delta-V and propellant budget for a continuous low-thrust
%   transfer between two circular orbits around a central body.  Combines:
%     - Edelbaum (1961) analytical delta-V — optimal continuous-thrust spiral
%     - Tangential-thrust RK4 propagation  — actual spiral trajectory
%     - Tsiolkovsky propellant accounting  — mass budget for given thruster
%
%   For delta_i = 0 the Edelbaum formula reduces to |v0 - v1|, the correct
%   result for a slow (many-revolution) circular spiral.  For combined orbit-
%   raising + plane change it gives the minimum-delta-V solution, which is
%   lower than the impulsive two-burn equivalent when delta_i is large.
%
%   The RK4 propagator uses prograde thrust (along velocity) for outward
%   spirals and retrograde thrust for inward spirals.  It is a good match
%   to Edelbaum for low thrust-to-weight ratios; at higher accelerations
%   the spiral departs from Edelbaum assumptions and the two will differ.
%
%   Inputs:
%     body      - central body struct from constants() (e.g., bodies.Earth)
%     r0_km     - initial circular orbit radius (km from body centre)
%     r1_km     - target  circular orbit radius (km from body centre)
%     options   (optional struct):
%       .inclinationChangeDeg  inclination change (deg)                [0]
%       .thrustN               engine thrust (N); 0 => Edelbaum only   [0]
%       .isp                   specific impulse (s)                   [3000]
%       .wetMass               spacecraft wet mass at start (kg)      [1000]
%       .nStepsPerOrbit        RK4 steps per osculating orbital period  [50]
%       .nOutputPoints         decimated trajectory output points     [2000]
%
%   Outputs (result struct):
%     .deltaV            Edelbaum delta-V (km/s)
%     .tof               estimated time of flight (s)  [NaN if thrustN = 0]
%     .tofDays           same in days
%     .propellantMass    propellant consumed (kg)        [NaN if thrustN = 0]
%     .finalMass         final spacecraft mass (kg)      [NaN if thrustN = 0]
%     .details           struct with component values
%     .trajectory        struct with decimated RK4 history (empty if thrustN = 0)
%                          .t  .x  .y  .vx  .vy  .r  .speed  .mass

if nargin < 4, options = struct(); end
if ~isfield(options,'inclinationChangeDeg'), options.inclinationChangeDeg = 0;    end
if ~isfield(options,'thrustN'),              options.thrustN              = 0;    end
if ~isfield(options,'isp'),                  options.isp                  = 3000; end
if ~isfield(options,'wetMass'),              options.wetMass              = 1000; end
if ~isfield(options,'nStepsPerOrbit'),       options.nStepsPerOrbit       = 50;   end
if ~isfield(options,'nOutputPoints'),        options.nOutputPoints        = 2000; end

mu  = body.mu;        % km^3/s^2
g0  = 9.80665e-3;     % km/s^2
F   = options.thrustN * 1e-3;   % N -> kN  (km unit system: 1 kN = 1 kg*km/s^2)
Isp = options.isp;
m0  = options.wetMass;

%% ---- Edelbaum analytical delta-V ------------------------------------
v0 = sqrt(mu / r0_km);   % initial circular speed (km/s)
v1 = sqrt(mu / r1_km);   % final   circular speed (km/s)

deltaI_rad = options.inclinationChangeDeg * pi / 180;
% Edelbaum (1961): dV = sqrt(v0^2 + v1^2 - 2*v0*v1*cos(pi*deltaI/2))
%   delta_i in radians; for delta_i = 0 reduces to |v0 - v1|
deltaV = sqrt(v0^2 + v1^2 - 2*v0*v1*cos(pi * deltaI_rad / 2));

%% ---- Tsiolkovsky propellant accounting -------------------------------
if F > 0
    mf             = m0 * exp(-deltaV / (Isp * g0));
    propellantMass = m0 - mf;
    tof            = propellantMass * (Isp * g0) / F;   % s (constant thrust)
else
    mf             = NaN;
    propellantMass = NaN;
    tof            = NaN;
end

%% ---- Pack result -----------------------------------------------------
result               = struct();
result.deltaV        = deltaV;
result.tof           = tof;
result.tofDays       = tof / 86400;
result.propellantMass = propellantMass;
result.finalMass     = mf;
result.details       = struct( ...
    'v0',                   v0,       ...
    'v1',                   v1,       ...
    'r0',                   r0_km,    ...
    'r1',                   r1_km,    ...
    'inclinationChangeDeg', options.inclinationChangeDeg, ...
    'thrust',               options.thrustN, ...
    'isp',                  Isp,      ...
    'wetMass',              m0);
result.trajectory    = struct();

if F <= 0
    return
end

%% ---- Tangential-thrust RK4 spiral -----------------------------------
% State vector (2-D Cartesian): [x; y; vx; vy; mass]
% Thrust direction: prograde (+F) for outward spiral, retrograde (-F) inward.
thrustDir = sign(r1_km - r0_km);
if thrustDir == 0, thrustDir = 1; end

state = [r0_km; 0; 0; v0; m0];

% Estimate step-decimation ratio.
% nStepsEst uses the initial (smallest) step size, so it overestimates
% the true step count for outward spirals. Cap stepsPerOut so we always
% retain at least minPtsPerOrbit stored points per orbit — otherwise the
% spiral looks like a polygon.
T0             = 2*pi * sqrt(r0_km^3 / mu);
dt0            = T0 / options.nStepsPerOrbit;
nStepsEst      = max(1, ceil(tof / dt0));
minPtsPerOrbit = 20;
maxStepsPerOut = max(1, floor(options.nStepsPerOrbit / minPtsPerOrbit));
stepsPerOut    = min(maxStepsPerOut, max(1, floor(nStepsEst / options.nOutputPoints)));

% Pre-allocate buffer sized for the actual expected output
nMax = max(options.nOutputPoints, ceil(nStepsEst / stepsPerOut)) + 500;
buf  = zeros(7, nMax);   % rows: [t; x; y; vx; vy; mass; r]

nOut       = 1;
buf(:, 1)  = [0; state; r0_km];
t          = 0;
stepCount  = 0;

while true
    r_cur = sqrt(state(1)^2 + state(2)^2);

    % Adaptive step: fraction of current osculating orbit period
    T_cur = 2*pi * sqrt(r_cur^3 / mu);
    dt    = T_cur / options.nStepsPerOrbit;
    dt    = min(dt, max(0, tof - t) + 1e-6);   % don't overshoot budget

    % RK4
    k1 = dt * ltODE(state, mu, F * thrustDir, Isp, g0);
    k2 = dt * ltODE(state + 0.5*k1, mu, F * thrustDir, Isp, g0);
    k3 = dt * ltODE(state + 0.5*k2, mu, F * thrustDir, Isp, g0);
    k4 = dt * ltODE(state +     k3, mu, F * thrustDir, Isp, g0);
    state = state + (k1 + 2*k2 + 2*k3 + k4) / 6;

    % Clamp mass to final mass floor
    if state(5) < mf, state(5) = mf; end

    t         = t + dt;
    stepCount = stepCount + 1;
    r_new     = sqrt(state(1)^2 + state(2)^2);

    % Store decimated point
    if mod(stepCount, stepsPerOut) == 0
        if nOut < nMax
            nOut = nOut + 1;
            buf(:, nOut) = [t; state; r_new];
        end
    end

    % Termination: reached target radius, propellant exhausted, or over budget
    outward = (thrustDir > 0);
    reachedTarget = (outward && r_new >= r1_km) || (~outward && r_new <= r1_km);
    done = reachedTarget || state(5) <= mf * 1.00001 || t >= tof * 1.002;
    if done
        if nOut < nMax
            nOut = nOut + 1;
            buf(:, nOut) = [t; state; r_new];
        end
        break
    end
end

buf = buf(:, 1:nOut);
result.trajectory = struct( ...
    't',     buf(1,:),                           ...
    'x',     buf(2,:),                           ...
    'y',     buf(3,:),                           ...
    'vx',    buf(4,:),                           ...
    'vy',    buf(5,:),                           ...
    'mass',  buf(6,:),                           ...
    'r',     buf(7,:),                           ...
    'speed', sqrt(buf(4,:).^2 + buf(5,:).^2));
end


%% ---- ODE: tangential low-thrust in 2-D Cartesian ---------------------
function ds = ltODE(state, mu, F_kN, Isp, g0)
% F_kN positive => prograde (along velocity), negative => retrograde

x  = state(1);  y  = state(2);
vx = state(3);  vy = state(4);
m  = max(state(5), 1e-9);

r3 = (x^2 + y^2)^1.5;
v  = sqrt(vx^2 + vy^2);

% Gravitational acceleration
ax_g = -mu * x / r3;
ay_g = -mu * y / r3;

% Thrust along (or against) velocity direction
if v > 1e-10
    ax_t = F_kN * vx / (m * v);
    ay_t = F_kN * vy / (m * v);
    mdot = -abs(F_kN) / (Isp * g0);
else
    ax_t = 0;  ay_t = 0;  mdot = 0;
end

ds = [vx; vy; ax_g + ax_t; ay_g + ay_t; mdot];
end
