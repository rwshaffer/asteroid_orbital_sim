function [delta_v, TOF] = escape_to_asteroid_plane(target,earth,num_phasing_orbits,parking_orbit,start_date)
% Ryan Shaffer, AU Diggers
% Interplanetary trajectory from Earth to a target body.
% Uses the patched conic approximation.
% Neglects target gravity (assuming asteroid target with small mass).
% Currently only configured for three-burn impulsive solutions (see trajectory steps below).
%
% INPUTS:
% target_orbit: struct containing Keplerian elements of target around Sun
%   target.epoch - epoch at which orbit is based (JD)
%   target.sma - semimajor axis (km)
%   target.ecc - eccentricity
%   target.inc - inclination (rad)
%   target.raan - right ascension of the ascending node (rad)
%   target.aop - argument of periapsis (rad)
%   target.nu - true anomaly (rad)
% earth: struct containing Keplerian elements of Earth around Sun
%   contains the same fields as target struct above
% parking_orbit: initial satellite insertion: "GTO" or "LEO". Default GTO.
% num_phasing_orbits: Number of heliocentric orbits used to match target phase.
% start_date: datetime variable - date to begin looking for launch windows. Default 1/1/2023.
%   use the following function to create a datetime variable: datetime(year,month,day)
%   where year, month, day are integers
%
% OUTPUTS:
% delta_v - Total delta-V over the full trajectory
% TOF - Total time of flight from start to finish
%
%
% TRAJECTORY PROFILE:
% Step 0: Maneuver from whatever launch orbit to parking orbit that is coplanar to asteroid orbit
%   (i.e., the orbit is inclined to match the asteroid's orbit)
% Step 1: At target's ascending or descending node, burn to Earth Escape
% Step 2: Heliocentric Hohmann transfer to the other node
% Step 3: Phasing maneuver to intercept target after 1 to 2 orbits
% Step 4: Match target velocity and rendezvous


%% HANDLE INPUTS
if nargin < 5
    % start_date not given - default is 1/1/2023
    start_date = datetime(2023,1,1);
    if nargin < 4
        % parking orbit not given - default is GTO
        parking_orbit = "GTO";
        if nargin < 3
            % num_phasing_orbits not given - default 1
            num_phasing_orbits = 1;
        end
    end
end

% Convert start_date from datetime to Julian date
start_date = juliandate(start_date);

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

%% ORBIT SETUP

% Starting satellite orbit parameters (around Earth)
switch parking_orbit
    case "GTO"
        parking.rp = earth_rad + 400; % Radius of perigee for parking orbit (km)
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

%% ASCENDING/DESCENDING NODES
% Satellite will escape from Earth when the Earth passes either the target's ascending or descending node.

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


%% STEP 1a: Find launch window when Earth first crosses target's line of nodes
parking.desc_wait_time = TimeOfFlight(earth,sun_mu,earth.nu_target_desc); % time (sec) spent waiting to reach true anomaly of desc. node
parking.asc_wait_time = TimeOfFlight(earth,sun_mu,earth.nu_target_asc); % time (sec) spent waiting to reach true anomaly of ascending node
if parking.desc_wait_time < parking.asc_wait_time
    % First time Earth crosses the line of nodes is at the descending node
    launch_node = 'descending';
    parking.wait_time = parking.desc_wait_time; % time (sec) spent waiting to reach Earth escape window
    earth.nu_launch = earth.nu_target_desc; % Earth true anomaly at launch position
    earth.rad_launch = earth.rad_desc; % Earth orbital radius at launch position
    target.rad_launch = target.rad_desc; % Target's orbital radius at node near launch position
    target.rad_opp_launch = target.rad_asc; % Target's orbital radius at node opposite of launch position
else
    % First time Earth crosses the line of nodes is at the ascending node
    launch_node = 'ascending';
    parking.wait_time = parking.asc_wait_time; % time (sec) spent waiting to reach Earth escape window
    earth.nu_launch = earth.nu_target_asc; % Earth true anomaly at launch position
    earth.rad_launch = earth.rad_asc; % Earth orbital radius at launch position
    target.rad_launch = target.rad_asc; % Target's orbital radius at node near launch position
    target.rad_opp_launch = target.rad_desc; % Target's orbital radius at node opposite of launch position
end

launch_window = earth.epoch + parking.wait_time * sec2day; % Julian day of launch window (node crossing)
earth.nu = earth.nu_launch; % Update Earth true anomaly to launch position
earth.epoch = launch_window; % Update Earth epoch to launch window

target.wait_TOF = (launch_window - target.epoch) * day2sec; % Time of flight (sec) between target epoch and launch window 
target.nu_launch = KeplerProblem(target,sun_mu,target.wait_TOF); % Target true anomaly at launch time (rad)
target.nu = target.nu_launch; % Update true anomaly to position at launch
target.epoch = launch_window; % Update target epoch to launch window

%% Step 2: Heliocentric Hohmann transfer 
%   from Earth orbit at intersection of target line of nodes to target's other node
transfer.sma = 1/2 * (earth.rad_launch + target.rad_opp_launch); % transfer orbit semi-major axis (km)
transfer.eps = -sun_mu/(2*transfer.sma); % transfer orbit energy (km^2/s^2)
transfer.v1_rel_sun = sqrt(2*(transfer.eps + sun_mu/earth.rad_launch)); % velocity along transfer orbit at edge of Earth's SOI (km/s)
transfer.v2 = sqrt(2*(transfer.eps + sun_mu/target.rad_opp_launch)); % velocity along transfer orbit at helio. apoapsis (km/s)
earth.v1 = sqrt(2*(earth.eps + sun_mu/earth.rad_launch)); % Earth's velocity at point of satellite departure (km/s)
transfer.phi_1 = 0; % satellite asymptote angle (rad)
v1_rel_earth = sqrt(transfer.v1_rel_sun^2 + earth.v1^2 - 2*transfer.v1_rel_sun * earth.v1 * cos(transfer.phi_1)); % satellite velocity relative to Earth at edge of Earth's SOI (km/s)

%% Step 1b: Hyperbolic escape from Earth's SOI
escape.v_inf = v1_rel_earth; % Hyperbolic escape velocity from Earth = satellite relative velocity at SOI
escape.C3 = escape.v_inf^2; % Hyperbolic characteristic energy (km^2/s^2)
delta_v.escape = C3_to_dV(escape.C3,parking_orbit); % delta-V required to escape Earth's SOI as desired (km/s)

%% Step 3: Phasing burn at other node to meet target
% Step 3a: Find phase angle between spacecraft and target
transfer.TP = sqrt(4*pi^2/sun_mu * transfer.sma^3); % Orbital period (sec) for Heliocentric Hohmann transfer (step 2)
transfer.TOF = transfer.TP/2; % Time of Flight (sec) for spacecraft to complete Hohmann transfer half-orbit
target.nu = KeplerProblem(target,sun_mu,transfer.TOF); % target true anomaly (rad) when S/C reaches node
target.epoch = target.epoch + transfer.TOF * sec2day; % Update target epoch to match new true anomaly
delta_phi = target.nu_asc - target.nu; % Phase angle (rad) between S/C and target when S/C reaches the target's node

% Step 3b: Time of Flight for target to travel that phase angle
phasing.time = TimeOfFlight(target,sun_mu,target.nu_asc); % time (sec) of phasing required

% Step 3c: Phasing orbit needed for S/C to rendezvous after K orbits
K = num_phasing_orbits; % Number of complete orbits in phasing trajectory.
phasing.TP = target.TP + phasing.time/K; % Orbital period (sec) required for phasing orbit
phasing.TOF = phasing.TP * K; % Time of flight for phasing orbit (sec)
phasing.sma = (phasing.TP^2 * sun_mu/(4*pi^2))^(1/3); % semi-major axis (km) for phasing orbit
if phasing.sma > target.rad_opp_launch
    phasing.rp = target.rad_opp_launch; % node opposite launch (location where burn is made) is periapsis of phase orbit
    phasing.ra = 2*phasing.sma - phasing.rp; % radius of apoapsis (km)
else
    phasing.ra = target.rad_opp_launch; % node opposite launch (location where burn is made) is apoapsis of phase orbit
    phasing.rp = 2*phasing.sma - phasing.ra; % radius of periapsis (km)
end
phasing.ecc = (phasing.ra - phasing.rp)/(phasing.ra + phasing.rp); % eccentricity of phase orbit
phasing.p = phasing.sma * (1 - phasing.ecc^2); % phasing orbit parameter
phasing.h = sqrt(phasing.p * sun_mu); % phasing orbit angular momentum
if phasing.sma > target.rad_opp_launch
    phasing.vp = phasing.h/phasing.rp; % velocity at periapsis (burn) on phase orbit
    phasing.v3 = phasing.vp; % S/C velocity (km/s) immediately after phasing burn
else
    phasing.va = phasing.h/phasing.ra; % velocity at apoapsis (burn) on phase orbit
    phasing.v3 = phasing.va; % S/C velocity (km/s) immediately after phasing burn
end
delta_v.phasing = abs(phasing.v3 - transfer.v2); % delta-V (km/s) required to enter correct phasing orbit (km/s)

%% Step 4: Match target velocity to rendezvous
arrival_epoch = target.epoch + phasing.TOF * sec2day;
phasing.v4 = phasing.v3; % point of rendezvous is still target's ascending node, so velocity is unchanged
target.v4 = sqrt(2*(target.eps + sun_mu/target.rad_opp_launch)); % target's velocity at rendezvous (node opposite launch) (km/s)
delta_v.rdvs = abs(phasing.v4 - target.v4); % delta-V (km/s) required to match target velocity

%% SUMMARY
% Total Delta-V
delta_v.total = delta_v.escape + delta_v.phasing + delta_v.rdvs;
fprintf("Burn 1 (escape): \t %0.3f km/s \n", delta_v.escape)
fprintf("Burn 2 (phasing): \t %0.3f km/s \n", delta_v.phasing)
fprintf("Burn 3 (rendezvous): %0.3f km/s \n", delta_v.rdvs)
fprintf("Total Delta-V: \t %0.3f km/s \n",delta_v.total)
% Total time of flight
TOF = transfer.TOF + phasing.TOF; % Total mission time of flight (sec). Doesn't include Earth escape or wait times - include either?
fprintf("Time of Flight: %0.2f days / %0.3f years \n",TOF*sec2day, TOF*sec2day/365)

% Launch date (from parking orbit)
fprintf("Launch Epoch: %s \n",string(datetime(launch_window,'convertfrom','juliandate')))
% Arrival date to target
fprintf("Arrival Epoch: %s \n",string(datetime(arrival_epoch,'convertfrom','juliandate')))
fprintf("\nNote: still need to incorporate maneuvers around Earth to get into inclined parking orbit. \n\n")



%% Organize trajectory values (delta V's, important times in trajectory, orbital elements, etc. 
% for use in STK simulation









end




