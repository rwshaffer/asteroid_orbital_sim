function [delta_v, TOF] = interplanetary_impulsive(target,parking_orbit,start_date,earth_orbit)
% Ryan Shaffer, AU Diggers
% Interplanetary trajectory from Earth to a target body.
% Uses the patched conic approximation.
% Neglects target gravity (assuming asteroid target with small mass).
% Currently only configured for three-burn impulsive solutions (see trajectory steps below).
%
% INPUTS:
% target_orbit: struct containing Keplerian elements of target aronud sun
%   target.epoch - epoch at which orbit is based (JD)
%   target.sma - semimajor axis (km)
%   target.ecc - eccentricity
%   target.inc - inclination (rad)
%   target.raan - right ascension of the ascending node (rad)
%   target.aop - argument of periapsis (rad)
%   target.nu - true anomaly (rad)
% parking_orbit: initial satellite insertion: "GTO" or "LEO". Default GTO.
% start_date: Julian day to begin looking for trajectories. Default 1/1/2022.
%   use the following Julian dates for future years:
%   1/1/2022: 2459580.5
%   1/1/2023: 2459945.5
%   1/1/2024: 2460310.5
%   1/1/2025: 2460676.5
%   https://www.aavso.org/jd-calculator - Online conversion between calendar and Julian dates
% earth_orbit: "circular" or "eccentric". Default to approximation of circular orbit.
%   note: eccentric earth_orbit not yet supported and probs won't be. sucks

if nargin < 4
    earth_orbit = "circular";
    if nargin < 3
        start_date = 2459580.5;
        if nargin < 2
            parking_orbit = "GTO";
        end
    end
end

%% UNIT CONVERSIONS
au2km = 1.496e8; % Astronomical units to kilometers
deg2rad = pi/180; % Degrees to radians
rad2deg = 180/pi; % Radians to degrees
day2sec = 86400; % Days to seconds
sec2day = 1/86400; % Seconds to days

%% CONSTANTS
sun_mu = 1.3271233e11; % Gravitational parameter (km^3/s^2)
earth_mu = 3.986e5; % Gravitational parameter (km^3/s^2)
earth_rad = 6378; % Radius of Earth (km)

% Earth's orbit around the Sun:
earth.epoch = 2451545; % epoch at which orbit is based (JD)
earth.sma = 1.00000011 * au2km; % semi-major axis (km)
if earth_orbit == "eccentric"
    earth.ecc = 0.01671; % eccentricity
    earth.inc = 0 * deg2rad; % inclination above ecliptic (rad)
    earth.raan = 197.91 * deg2rad; % right ascension (rad)
    earth.aop = 267.452 * deg2rad; % argument of periapsis (rad)
    j2000_mean_long = 100.46435 * deg2rad; % Earth's mean longitude at epoch (rad)
    mean_anomaly = mod(j2000_mean_long - (earth.raan + earth.aop),2*pi); % mean anomaly (rad)
    earth.nu = convert_anomalies(mean_anomaly,earth.ecc,'mean'); % true anomaly (rad)
else
    earth.ecc = 0; % eccentricity
    earth.inc = 0; % inclination above ecliptic (rad)
    earth.raan = 0; % right ascension (rad)
    earth.aop = 0; % argument of periapsis (rad)
    j2000_mean_long = 100.46435 * deg2rad; % Earth's mean longitude at epoch (rad)
    earth.nu = j2000_mean_long; % true anomaly (rad) w.r.t vernal equinox
end

% Starting satellite orbit parameters (around Earth)
switch parking_orbit
    case "GTO"
        parking.rp = earth_rad + 200; % Radius of perigee for parking orbit (km)
        parking.ra = 42164; % Radius of apogee for parking orbit (km)
    case "LEO"
        parking.rp = earth_rad + 200; % Radius of perigee for parking orbit (km)
        parking.ra = parking.rp; % Radius of apogee for parking orbit (km)
end
parking.sma = (parking.rp + parking.ra)/2; % semi-major axis (km)
parking.ecc = (parking.ra - parking.rp)/(parking.ra + parking.rp); % eccentricity
parking.eps = -earth_mu/(2*parking.sma); % parking orbital energy (km^2/s^2)


% Earth orbital parameters
earth.eps = -sun_mu/(2*earth.sma); % Earth orbital energy (km^2/s^2)
earth.p = earth.sma * (1 - earth.ecc^2); % Parameter (used in trajectory equation)
earth.traj = @(nu) earth.p/(1 + earth.ecc*cos(nu)); % Trajectory equation, input true anomaly (rad). Output rad (km).

% Target object orbital parameters
target.eps = -sun_mu/(2*target.sma); % target orbital energy (km^2/s^2)
target.TP = sqrt(4*pi^2/sun_mu * target.sma^3); % Orbital period (sec) for target
target.p = target.sma * (1 - target.ecc^2); % Parameter (used in trajectory equation)
target.traj = @(nu) target.p/(1 + target.ecc*cos(nu)); % Trajectory equation, input true anomaly (rad). Output rad (km).

% Propagate Earth and target to same starting epoch (start_date)
earth_TOF = (start_date - earth.epoch) * day2sec; % Time (sec) between Earth's epoch and desired epoch
earth.nu = KeplerProblem(earth,sun_mu,earth_TOF); % True anomaly (rad) at new epoch
earth.epoch = start_date; % Update Earth epoch
target_TOF = (start_date - target.epoch) * day2sec; % Time (sec) between target's epoch and desired epoch
target.nu = KeplerProblem(target,sun_mu,target_TOF); % True anomaly (rad) at new epoch
target.epoch = start_date; % Update target epoch

%% GENERAL HELIOCENTRIC TRAJECTORY SETUP
% Trajectory steps:
% Step 0: Maneuver from whatever launch orbit to parking orbit that is coplanar to asteroid orbit
%   (i.e., the orbit is inclined to match the asteroid's orbit)
% Step 1: At target's ascending or descending node, burn to Earth Escape
% Step 2: Heliocentric Hohmann transfer to the other node
% Step 3: Phasing maneuver to intercept target after 1 to 2 orbits
% Step 4: Match target velocity and rendezvous

% Find ascending and descending nodes for target orbit
target.nu_asc = 2*pi - target.aop; % true anomaly of ascending node (rad)
target.nu_desc = target.nu_asc + pi; % true anomaly of descending node (rad)
earth.nu_target_asc = target.raan; % Earth's true anomaly at target ascending node (rad)
earth.nu_target_desc = earth.nu_target_asc + pi; % Earth's true anomaly at target descending node (rad)

% Radii of target orbit at each node
target.rad_asc = target.traj(target.nu_asc);
target.rad_desc = target.traj(target.nu_desc);
% Radius of Earth orbit at each of the target's nodes
earth.rad_asc = earth.traj(earth.nu_target_asc);
earth.rad_desc = earth.traj(earth.nu_target_desc);


%% STEP 1a: Find launch window when Earth crosses target's line of nodes
parking.wait_time = TimeOfFlight(earth,sun_mu,earth.nu_target_desc); % time (sec) spent waiting to reach true anomaly
launch_window = earth.epoch + parking.wait_time * sec2day; % Julian day of launch window (node crossing)
earth.nu = earth.nu_target_desc; % Update Earth true anomaly to launch position
earth.epoch = launch_window; % Update Earth epoch to launch window

target.wait_TOF = (launch_window - target.epoch) * day2sec; % Time of flight (sec) between target epoch and launch window 
target.nu_launch = KeplerProblem(target,sun_mu,target.wait_TOF); % Target true anomaly at launch time (rad)
target.nu = target.nu_launch; % Update true anomaly to position at launch
target.epoch = launch_window; % Update target epoch to launch window

%% Step 2: Heliocentric Hohmann transfer 
%   from Earth orbit (near target desc. node) to target ascending node
transfer.sma = 1/2 * (earth.rad_desc + target.rad_asc); % transfer orbit semi-major axis (km)
transfer.eps = -sun_mu/(2*transfer.sma); % transfer orbit energy (km^2/s^2)
transfer.v1_rel_sun = sqrt(2*(transfer.eps + sun_mu/earth.rad_desc)); % velocity along transfer orbit at edge of Earth's SOI (km/s)
transfer.v2 = sqrt(2*(transfer.eps + sun_mu/target.rad_asc)); % velocity along transfer orbit at helio. apoapsis (km/s)
earth.v1 = sqrt(2*(earth.eps + sun_mu/earth.rad_desc)); % Earth's velocity at point of satellite departure (km/s)
transfer.phi_1 = 0; % satellite asymptote angle (rad)
v1_rel_earth = sqrt(transfer.v1_rel_sun^2 + earth.v1^2 - 2*transfer.v1_rel_sun * earth.v1 * cos(transfer.phi_1)); % satellite velocity relative to Earth at edge of Earth's SOI (km/s)

%% Step 1b: Hyperbolic escape from Earth's SOI
escape.v_inf = v1_rel_earth; % Hyperbolic escape velocity from Earth = satellite relative velocity at SOI
escape.eps = escape.v_inf^2/2; % Hyperbolic trajectory energy (km^2/s^2)
escape.r0 = parking.rp; % Departure radius (km). Assume departure at perigee
escape.v0 = sqrt(2*(escape.eps + earth_mu/escape.r0)); % Velocity of sat. on transfer orbit at departure point
parking.v0 = sqrt(2*(parking.eps + earth_mu/escape.r0)); % Velocity of sat. in parking orbit at departure point
delta_v1 = abs(escape.v0 - parking.v0); % delta-V required to escape Earth's SOI as desired (km/s)

%% Step 3: Phasing burn at other node to meet target
% Step 3a: Find phase angle between spacecraft and target
transfer.TP = sqrt(4*pi^2/sun_mu * transfer.sma^3); % Orbital period (sec) for Heliocentric Hohmann transfer (step 2)
transfer.TOF = transfer.TP/2; % Time of Flight (sec) for spacecraft to complete half of Hohmann transfer orbit
target.nu = KeplerProblem(target,sun_mu,transfer.TOF); % target true anomaly (rad) when S/C reaches node
target.epoch = target.epoch + transfer.TOF * sec2day; % Update target epoch to match new true anomaly
delta_phi = target.nu_asc - target.nu; % Phase angle (rad) between S/C and target when S/C reaches the target's node

% Step 3b: Time of Flight for target to travel that phase angle
phasing.time = TimeOfFlight(target,sun_mu,target.nu_asc); % time (sec) of phasing required

% Step 3c: Phasing orbit needed for S/C to rendezvous after K orbits
K = 2; % Number of complete orbits in phasing trajectory. (Make this an input?)
phasing.TP = target.TP + phasing.time/K; % Orbital period (sec) required for phasing orbit
phasing.TOF = phasing.TP * K; % Time of flight for phasing orbit (sec)
phasing.sma = (phasing.TP^2 * sun_mu/(4*pi^2))^(1/3); % semi-major axis (km) for phasing orbit
if phasing.sma > target.rad_asc
    phasing.rp = target.rad_asc; % ascending node (location where burn is made) is periapsis of phase orbit
    phasing.ra = 2*phasing.sma - phasing.rp; % radius of apoapsis (km)
else
    phasing.ra = target.rad_asc; % ascending node (location where burn is made) is apoapsis of phase orbit
    phasing.rp = 2*phasing.sma - phasing.ra; % radius of periapsis (km)
end
phasing.ecc = (phasing.ra - phasing.rp)/(phasing.ra + phasing.rp); % eccentricity of phase orbit
phasing.p = phasing.sma * (1 - phasing.ecc^2); % phasing orbit parameter
phasing.h = sqrt(phasing.p * sun_mu); % phasing orbit angular momentum
if phasing.sma > target.rad_asc
    phasing.vp = phasing.h/phasing.rp; % velocity at periapsis (burn) on phase orbit
    phasing.v3 = phasing.vp; % S/C velocity (km/s) immediately after phasing burn
else
    phasing.va = phasing.h/phasing.ra; % velocity at apoapsis (burn) on phase orbit
    phasing.v3 = phasing.va; % S/C velocity (km/s) immediately after phasing burn
end
delta_v2 = abs(phasing.v3 - transfer.v2); % delta-V (km/s) required to enter correct phasing orbit (km/s)

%% Step 4: Match target velocity to rendezvous
arrival_epoch = target.epoch + phasing.TOF * sec2day;
phasing.v4 = phasing.v3; % point of rendezvous is still target's ascending node, so velocity is unchanged
target.v4 = sqrt(2*(target.eps + sun_mu/target.rad_asc)); % target's velocity at rendezvous (ascending node) (km/s)
delta_v3 = abs(phasing.v4 - target.v4); % delta-V (km/s) required to match target velocity

%% SUMMARY
% Total Delta-V
delta_v = delta_v1 + delta_v2 + delta_v3;
fprintf("Burn 1 (escape): \t %0.3f km/s \n", delta_v1)
fprintf("Burn 2 (phasing): \t %0.3f km/s \n", delta_v2)
fprintf("Burn 3 (rendezvous): %0.3f km/s \n", delta_v3)
fprintf("Total Delta-V: \t %0.3f km/s \n",delta_v)
% Total time of flight
TOF = transfer.TOF + phasing.TOF; % Total mission time of flight (sec). Doesn't include Earth escape or wait times - include either?
fprintf("Time of Flight: %0.2f days / %0.3f years \n",TOF*sec2day, TOF*sec2day/365)

% Launch date (from parking orbit)
fprintf("Launch Epoch: %0.0f \n",launch_window)
% Arrival date to target
fprintf("Arrival Epoch: %0.0f \n",arrival_epoch)
fprintf("\nNote: still need to incorporate maneuvers around Earth to get into inclined parking orbit. \n\n")
end
