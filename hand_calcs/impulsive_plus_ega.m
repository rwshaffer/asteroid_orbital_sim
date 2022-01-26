function [delta_v, TOF] = impulsive_plus_ega(target,parking_orbit,num_phasing_orbits,start_date,earth_orbit)
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
% num_phasing_orbits: Number of heliocentric orbits used to match target phase.
% start_date: Julian day to begin looking for trajectories. Default 1/1/2022.
%   use the following Julian dates for future years:
%   1/1/2022: 2459580.5
%   1/1/2023: 2459945.5
%   1/1/2024: 2460310.5
%   1/1/2025: 2460676.5
%   https://www.aavso.org/jd-calculator - Online conversion between calendar and Julian dates
% earth_orbit: "circular" or "eccentric". Default to approximation of circular orbit.
%   note: eccentric earth_orbit not yet supported and probs won't be. sucks

if nargin < 5
    earth_orbit = "circular";
    if nargin < 4
        start_date = 2459580.5;
        if nargin < 3
            num_phasing_orbits = 1;
            if nargin < 2
                parking_orbit = "GTO";
            end
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
earth.TP = sqrt((4*pi^2/sun_mu)*earth.sma^3); % Orbital period (sec)

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
% Step 0: Maneuver from whatever launch orbit to parking orbit (tbd)
% Step 1: Near target's ascending node, burn to Earth Escape with orbital period of nearly 1 year
% Step 2: Heliocentric Hohmann transfer to the descending node
% Step 3: Earth gravity assist at desc. node to match target inclination and RAAN
% Step 4: Burn at periapsis until apoapsis touches target's orbit
% Step 5: Burn at apoapsis to match target's speed at that position
% Step 6: Phasing maneuver to intercept target after 1 to 2 orbits
% Step 7: Match target velocity and rendezvous
% 
% Steps are completed out of order, based on the flow of information.

%% SETUP - Comparing Earth and target orbits
% Find ascending and descending nodes for target orbit
target.nu_asc = 2*pi - target.aop; % true anomaly of ascending node (rad)
target.nu_desc = target.nu_asc + pi; % true anomaly of descending node (rad)
earth.nu_target_asc = target.raan; % Earth's true anomaly when aligned with target ascending node (rad)
earth.nu_target_desc = earth.nu_target_asc + pi; % Earth's true anomaly when aligned with target descending node (rad)

% Radii of target orbit at each node
target.rad_asc = target.traj(target.nu_asc);
target.rad_desc = target.traj(target.nu_desc);
% Radius of Earth orbit when aligned with each of the target's nodes
earth.rad_asc = earth.traj(earth.nu_target_asc);
earth.rad_desc = earth.traj(earth.nu_target_desc);

%% Start with the Heliocentric Hohmann transfer, as it is the most constrained maneuver
%% STEP 2: Heliocentric Hohmann transfer from Earth departure at target's asc. node 
%   to EGA at target's desc. node
% This transfer orbit is defined by the need to launch from Earth and rendezvous with Earth
%   (constraining the period) and the optimal inclination and RAAN that will match the 
%   post-EGA orbital plane to the target's plane.
transfer1.TP = 366 * day2sec; % Orbital period (sec) must be near Earth's to intercept for EGA
transfer1.TOF = transfer1.TP/2; % Time of flight (sec) for S/C along heliocentric transfer orbit
transfer1.sma = (sun_mu/(4*pi^2) * (transfer1.TP)^2)^(1/3); % Semi-major axis (km)
transfer1.rp = earth.rad_asc; % Assuming leaving at the ascending node
transfer1.ecc = 1 - transfer1.rp/transfer1.sma; % eccentricity
transfer1.inc = 0.8 * deg2rad; % initial guess at inclination (rad) (optimized in step 3b)
transfer1.raan = target.raan; % initial guess at RAAN (rad) (optimized in step 3b)
transfer1.aop = 0; % Argument of periapsis (rad). Zero, since Earth escape burn is at ascending node.

%% STEP 3: Earth Gravity Assist (EGA) when Earth is at target's descending node
%% Step 3a: Propagate Earth and S/C to correct positions at target's desc. node
earth.TOF = TimeOfFlight(earth,sun_mu,earth.nu_target_desc); % Time of flight from current epoch (start_date) to descending node
earth.epoch = earth.epoch + earth.TOF * sec2day; % Update epoch to when Earth crosses target's desc. node
earth.nu = earth.nu_target_desc; % Update Earth's true anomaly
transfer1.epoch = earth.epoch; % S/C epoch matches Earth's at EGA
transfer1.nu = pi; % EGA occurs at S/C apoapsis (maybe not exactly true? Close approx.)

%% Step 3b: Find optimal initial orbit that will match target's orbit after EGA
% Note that changing the Hohmann transfer inclination and RAAN will not affect the other 
%   orbital elements, so step 3a is still valid.
transfer1 = optimize_EGA(transfer1,earth,target); % Change inc. and RAAN to match target inc. and RAAN

%% Step 3c: Compute final orbit after EGA
[transfer2,B_dot,rp_EGA] = earth_grav_assist(transfer1, earth); % Earth Gravity Assist
transfer2.epoch = transfer1.epoch; % Currently neglects flyby TOF.

%% Now that step 2 is complete,move backwards in time to 
%   step 1 to determine launch/Earth escape required 
%% Step 1a: Find launch window when Earth intersects S/C's Hohmann ascending node
% Back up time and Earth's true anomaly to start_date
time_backup = (start_date - earth.epoch)*day2sec; % Time (sec) from current epoch (at EGA) to start date. Should be negative
earth.nu = KeplerProblem(earth,sun_mu,time_backup); % True anomaly (rad) from propagating Earth's orbit backwards
earth.epoch = start_date; % Update Earth epoch to start date (JD)

% Find launch window - when Earth intersects the S/C heliocentric Hohmann transfer asc. node
earth.nu_sc_asc = transfer1.raan; % True anomaly (rad) where Earth crosses asc. node in Step 2 traj.
parking.wait_time = TimeOfFlight(earth,sun_mu,earth.nu_sc_asc); % time (sec) spent waiting to reach true anomaly
launch_window = earth.epoch + parking.wait_time * sec2day; % Julian day of launch window (node crossing)
earth.nu = earth.nu_sc_asc; % Update Earth true anomaly (rad) to launch position
earth.epoch = launch_window; % Update Earth epoch to launch window

% Propagate target to launch window position
target.wait_TOF = (launch_window - target.epoch) * day2sec; % Time of flight (sec) between target epoch and launch window 
target.nu_launch = KeplerProblem(target,sun_mu,target.wait_TOF); % Target true anomaly at launch time (rad)
target.nu = target.nu_launch; % Update true anomaly to position at launch
target.epoch = launch_window; % Update target epoch to launch window

%% Step 1b: Determine S/C velocity rel. to Earth at Earth escape
% Update S/C true anomaly to its value at launch
transfer1.nu_EGA = transfer1.nu; % Save true anomaly at EGA (may be needed later?)
transfer1.nu = 0; % True anomaly (rad) at Earth-escape, occurring at S/C periapsis.
% Obtain S/C 3D velocity components in sun-centered ecliptic frame from orbital elements:
[transfer1.r1_rel_sun,transfer1.v1_rel_sun] = Kepler2Carts(transfer1,'sun'); % Helio-Cartesian state of S/C at Earth escape
% Obtain Earth's 3D velocity components in sun-centered ecliptic frame from orbital elements:
[earth.r1_rel_sun,earth.v1_rel_sun] = Kepler2Carts(earth,'sun'); % Helio-Cartesian state of S/C at Earth escape
% S/C velocity relative to Earth in Earth-centered ecliptic (NOT EQUATORIAL) frame:
v1_rel_earth = transfer1.v1_rel_sun - earth.v1_rel_sun; % Velocity vector in Earth frame

%% Step 1c: Hyperbolic escape from Earth's SOI
% Currently a 2D solution, could be upgraded to 3D in order to work more closely with parking orbit
escape.v_inf = norm(v1_rel_earth); % Hyperbolic escape velocity from Earth = satellite relative velocity at SOI
escape.eps = escape.v_inf^2/2; % Hyperbolic trajectory energy (km^2/s^2)
escape.r0 = parking.rp; % Departure radius (km). Assume departure at perigee
escape.v0 = sqrt(2*(escape.eps + earth_mu/escape.r0)); % Velocity of sat. on transfer orbit at departure point
parking.v0 = sqrt(2*(parking.eps + earth_mu/escape.r0)); % Velocity of sat. in parking orbit at departure point
delta_v.escape = abs(escape.v0 - parking.v0); % Delta-V required to escape Earth's SOI as desired (km/s)

%% Jump back ahead in time to just after EGA. S/C is in target's orbital plane
%   and will now work on matching its orbit.
%% STEP 4: Burn at periapsis until apoapsis touches target's orbit
% Step 4a: Propagate to S/C periapsis
transfer2.nu_final = 0; % Transfer 2 orbit ends with a burn at periapsis (true anomaly = 0)
transfer2.TP = sqrt(4*pi^2/sun_mu * transfer2.sma^3); % Orbital period (sec) of transfer orbit
transfer2.TOF = TimeOfFlight(transfer2,sun_mu,transfer2.nu_final); % Time of flight (sec) from EGA to periapsis    
transfer2.epoch = transfer2.epoch + transfer2.TOF * sec2day; % Epoch (JD) at periapsis (when burn occurs)
transfer2.nu = transfer2.nu_final; % Update true anomaly to periapsis

% Step 4b: Find target orbital radius at S/C apoapsis
transfer2.loa = transfer2.raan + transfer2.aop + pi; % Longitude of apoapsis (rad) for S/C orbit
target.nu_SC_apo = transfer2.loa - target.raan - target.aop; % Target true anomaly (rad) aligned with S/C orbit apoapsis
target.rad_SC_apo = target.traj(target.nu_SC_apo); % Target radius when aligned with S/C apoapsis

% Step 4c: Calculate burn required to match radius of apoapsis with target radius at nu_SC_apo
transfer2.rp = transfer2.sma * (1 - transfer2.ecc); % S/C radius of periapsis (km)
transfer2.eps = -sun_mu/(2*transfer2.sma); % Pre-burn orbit (transfer2) specific energy (km^2/s^2)
transfer2.vp = sqrt(2*(transfer2.eps + sun_mu/transfer2.rp)); % Initial velocity at periapsis (km/s)

transfer3.rp = transfer2.rp; % Burn (transfer2 -> transfer3) occurs at periapsis
transfer3.ra = target.rad_SC_apo; % New S/C orbit (transfer3) meets target's orbit at apoapsis
transfer3.sma = (transfer3.rp + transfer3.ra)/2; % New semi-major axis (km)
transfer3.ecc = (transfer3.ra - transfer3.rp)/(transfer3.ra + transfer3.rp); % Eccentricity
transfer3.inc = transfer2.inc;
transfer3.raan = transfer2.raan;
transfer3.aop = transfer2.aop;
transfer3.epoch = transfer2.epoch; % Epoch at beginning of transfer 3 = epoch at end of transfer 2
transfer3.nu = transfer2.nu_final; % True anomaly (rad) does not change, since burn occurs at periapsis.
%    NOTE: this is not necessarily the case for asteroids orbiting inside Earth's orbit or ECA's. 
%    It is possible that the burn at periapsis would be retrograde enough to make it the new apoapsis.
%    Things like this must be fixed to make the code generalizable.
transfer3.eps = -sun_mu/(2*transfer3.sma); % Post-burn orbit (transfer3) specific energy (km^2/s^2)
transfer3.vp = sqrt(2*(transfer3.eps + sun_mu/transfer3.rp)); % Final velocity at periapsis (km/s)

delta_v.raise_apo = abs(transfer3.vp - transfer2.vp); % Delta-V required to raise the orbit's apoapsis to meet target orbit

%% STEP 5: Match target's orbital speed at transfer orbit apoapsis
% Once reaching its apoapsis, the S/C will be positioned along the target's orbit.
% A burn to match the target's velocity at that orbital position will place it in the same orbit.

% Step 5a: Propagate the S/C to apoapsis, and target to same epoch
transfer3.nu_final = pi; % Transfer 3 orbit ends with a burn at apoapsis (true anomaly = 0)
transfer3.TOF = TimeOfFlight(transfer3,sun_mu,transfer3.nu_final); % Time of flight (sec) from periapsis to apoapsis
transfer3.epoch = transfer3.epoch + transfer3.TOF * sec2day; % Epoch (JD) at periapsis (when burn occurs)
transfer3.nu = transfer3.nu_final; % Update true anomaly (rad) to final value

target.TOF = (transfer3.epoch - target.epoch) * day2sec; % Time of flight (sec) for target to catch up to current epoch

target.nu =  KeplerProblem(target,sun_mu,target.TOF); % Update target true anomaly (rad)
target.epoch = transfer3.epoch; % Update epoch to reflect current time

% Step 5b: S/C and target 3D velocities at apoapsis
[transfer3.pos_apo,transfer3.vel_apo] = Kepler2Carts(transfer3,'sun'); % Helio-Cartesian state of S/C at apoapsis
target.old_nu = target.nu; % Save current true anomaly (must change to an incorrect value next)
target.nu = target.nu_SC_apo; % (Briefly) change target's true anomaly to match the orbital position at S/C apoapsis
[target.pos_SC_apo,target.vel_SC_apo] = Kepler2Carts(target,'sun'); % Helio-Cartesian state of target at S/C apoapsis
target.nu = target.old_nu; % Replace original, correct true anomaly at current epoch

% Step 5c: Burn to match velocities
vel_SC_rel_target = target.vel_SC_apo - transfer3.vel_apo; % Velocity vector of target w.r.t. S/C
delta_v.orbit_match = norm(vel_SC_rel_target); % Delta-V required to match orbital velocities

% Step 5d: Update S/C Keplerian elements (should be same as target besides true anomaly)
transfer4 = Carts2Kepler(transfer3.pos_apo,target.vel_SC_apo,'sun'); % New orbital elements for S/C
transfer4.epoch = transfer3.epoch; % Epoch at start of transfer 4 = epoch at end of transfer 3
% Other orbital parameters of interest
transfer4.TP = sqrt(4*pi^2/sun_mu * transfer4.sma^3); % Orbital period (sec) of transfer orbit
transfer4.rp = transfer4.sma * (1 - transfer4.ecc); % Radius of periapsis (used in step 6)
transfer4.ra = transfer4.sma * (1 + transfer4.ecc); % Radius of apoapsis (used in step 6)
transfer4.eps = -sun_mu/(2*transfer4.sma); % Orbital specific energy (km^2/s^2)
transfer4.vp = sqrt(2*(transfer4.eps + sun_mu/transfer4.rp)); % Velocity of periapsis (km/s). Used in step 6.
transfer4.va = sqrt(2*(transfer4.eps + sun_mu/transfer4.ra)); % Velocity of periapsis (km/s). Used in step 6.

%% Step 6: Phasing burn at new apoapsis/periapsis (whichever is closer) to meet target
% Step 6a: Determine whether apoapsis or periapsis is coming first. Conduct phase burn at whichever is sooner.
if transfer4.nu > 180*deg2rad
    transfer4.burn_time = "periapsis";
    phasing.burn_rad = transfer4.rp; % Orbital radius at time of phasing burn
    transfer4.nu_final = 0; % Transfer4 will end (and phasing will begin) at periapsis (nu = 0)
    transfer4.v_final = transfer4.vp; % Velocity at instant before burn
else
    transfer4.burn_time = "apoapsis";
    phasing.burn_rad = transfer4.ra; % Orbital radius at time of phasing burn
    transfer4.nu_final = pi; % Transfer4 will end (and phasing will begin) at apoapsis (nu = 180 deg)
    transfer4.v_final = transfer4.va; % Velocity at instant before burn

end

% Step 6b: Propagate S/C and target to burn time (S/C reaching peri/apoapsis)
transfer4.TOF = TimeOfFlight(transfer4,sun_mu,transfer4.nu_final); % Time of flight (sec) from point of orbit matching to peri/apoapsis
transfer4.epoch = transfer4.epoch + transfer4.TOF * sec2day; % Epoch (JD) at time of phasing burn
transfer4.nu = transfer4.nu_final; % Update true anomaly

target.TOF = (transfer4.epoch - target.epoch) * day2sec; % Time of flight (sec) for target to catch up to current epoch
target.nu =  KeplerProblem(target,sun_mu,target.TOF); % Update target true anomaly (rad)
target.epoch = transfer4.epoch; % Update epoch to reflect current time

% Step 6c: Find phase angle between S/C and target, and TOF to travel that angle
phi = transfer4.nu - target.nu; % Phase angle (rad) between S/C and target
phi_travel = 2*pi - phi; % Angle (rad) traveled by target during one phasing orbit

phasing.phase_time = TimeOfFlight(target,sun_mu,target.nu + phi); % Time (sec) of phasing required

% Step 6d: Phasing orbit needed for S/C to rendezvous after K orbits
K = num_phasing_orbits; % Number of complete orbits in phasing trajectory (function input)
phasing.TP = target.TP + phasing.phase_time/K; % Orbital period (sec) required for phasing orbit
phasing.TOF = phasing.TP * K; % Time of flight for phasing orbit (sec)
phasing.sma = (phasing.TP^2 * sun_mu/(4*pi^2))^(1/3); % semi-major axis (km) for phasing orbit
phasing.inc = transfer4.inc; % Inclination unchanged
phasing.raan = transfer4.raan; % RAAN unchanged
phasing.aop = transfer4.aop; % arg. of periapsis unchanged (except in special cases in if statements below)
phasing.nu = transfer4.nu; % True anomaly unchanged (except in special cases in if statements below)
phasing.epoch = transfer4.epoch; % Update epoch (JD)
if phasing.sma > phasing.burn_rad
    phasing.rp = phasing.burn_rad; % burn location is periapsis of phase orbit
    phasing.ra = 2*phasing.sma - phasing.rp; % radius of apoapsis (km)
    phasing.nu_final = 0; % Phasing orbit ends at periapsis of phase orbit (where burn occurred)
    if transfer4.burn_time == "apoapsis" % Burn location has changed from apo to periapsis. Must update aop/nu
        phasing.nu = 0;
        phasing.aop = phasing.aop + pi;
    end
else
    phasing.ra = phasing.burn_rad; % burn location is apoapsis of phase orbit
    phasing.rp = 2*phasing.sma - phasing.ra; % radius of periapsis (km)
    phasing.nu_final = pi; % Phasing orbit ends at periapsis of phase orbit (where burn occurred)
    if transfer4.burn_time == "periapsis" % Burn location has changed from peri to apoapsis. Must update aop/nu
        phasing.nu = pi;
        phasing.aop = phasing.aop - pi;
    end
end
phasing.ecc = (phasing.ra - phasing.rp)/(phasing.ra + phasing.rp); % eccentricity of phase orbit
phasing.p = phasing.sma * (1 - phasing.ecc^2); % phasing orbit parameter
phasing.h = sqrt(phasing.p * sun_mu); % phasing orbit angular momentum
if phasing.sma > phasing.burn_rad
    phasing.vp = phasing.h/phasing.rp; % velocity at periapsis (burn) on phase orbit
    phasing.v3 = phasing.vp; % S/C velocity (km/s) immediately after phasing burn
else
    phasing.va = phasing.h/phasing.ra; % velocity at apoapsis (burn) on phase orbit
    phasing.v3 = phasing.va; % S/C velocity (km/s) immediately after phasing burn
end
delta_v.phasing = abs(phasing.v3 - transfer4.v_final); % Delta-V (km/s) required to enter correct phasing orbit

%% Step 7: Match target velocity to rendezvous
% Step 7a: Propagate S/C and target to rendezvous point
phasing.nu = phasing.nu_final; % True anomaly (rad) at end of phasing orbit, aka rendezvous
phasing.epoch = transfer4.epoch + phasing.TOF * sec2day; % Epoch (JD) at rendezvous
arrival_epoch = phasing.epoch; % Epoch at rendezvous

target.nu = KeplerProblem(target,sun_mu,phasing.TOF); % True anomaly (rad) of target at rendezvous (should match S/C)
target.epoch = phasing.epoch; % Update epoch to match rendezvous

% Step 7b: Cartesian states of S/C and target at rendezvous. Positions should match
[phasing.rdvs_pos,phasing.rdvs_vel] = Kepler2Carts(phasing,'sun'); % Cartesian state of S/C at end of phasing orbit
[target.rdvs_pos,target.rdvs_vel] = Kepler2Carts(target,'sun'); % Cartesian state of target at rendezvous

% Step 7c: Delta-V required to match velocities and rendezvous
vel_target_rel_SC = target.rdvs_vel - phasing.rdvs_vel; % Velocity vector of target rel. to S/C
delta_v.rdvs = norm(vel_target_rel_SC); % Delta-V (km/s) required to match velocities


%% SUMMARY
% Total Delta-V
delta_v.total = delta_v.escape + delta_v.raise_apo + delta_v.orbit_match + delta_v.phasing + delta_v.rdvs;
fprintf("Burn 1 (escape): \t %0.3f km/s \n", delta_v.escape)
fprintf("Burn 2 (raise apoapsis): \t %0.3f km/s \n", delta_v.raise_apo)
fprintf("Burn 3 (orbit matching): %0.3f km/s \n", delta_v.orbit_match)
fprintf("Burn 4 (phasing): %0.3f km/s \n", delta_v.phasing)
fprintf("Burn 5 (rendezvous): %0.3f km/s \n", delta_v.rdvs)
fprintf("Total Delta-V: \t %0.3f km/s \n",delta_v.total)
% Total time of flight
TOF = transfer1.TOF + transfer2.TOF + transfer3.TOF + transfer4.TOF + phasing.TOF; % Total mission time of flight (sec). 
% Doesn't include Earth escape time or EGA time of flight right now
fprintf("\nTime of Flight: %0.2f days / %0.3f years \n",TOF*sec2day, TOF*sec2day/365)

% Launch date (from parking orbit)
fprintf("Launch Epoch: %0.0f \n",launch_window)
% Arrival date to target
fprintf("Arrival Epoch: %0.0f \n\n",arrival_epoch)
%fprintf("\nNote: still need to incorporate maneuvers around Earth to get into inclined parking orbit. \n\n")
end
