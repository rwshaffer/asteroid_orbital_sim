function [delta_v, TOF] = ega_return(target,earth,start_date)
% Gen Gemond, AU Diggers
% Interplanetary trajectory from Earth to a target body.
% Uses the patched conic approximation.
% Neglects target gravity (assuming asteroid target with small mass).
% Currently only configured for impulsive solutions (see trajectory steps below).
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
% NOTE: in this case "target" denotes the asteroid (despite the satellite actually targeting Earth)
% earth: struct containing Keplerian elements of Earth around Sun
%   contains the same fields as target struct above
% start_date: datetime variable - date to begin looking for launch windows from asteroid. Default 1/1/2023.
%   use the following function to create a datetime variable: datetime(year,month,day)
%   where year, month, day are integers
%
% OUTPUTS:
% delta_v - Total delta-V over the full trajectory
% TOF - Total time of flight from start to finish
%
% TRAJECTORY PROFILE:
% Step 1: Wait until Earth-target phase angle is correct to perform a Hohmann transfer from asteroid to Earth
% Step 2: Perform Hohmann transfer and intercept Earth (but do not perform second burn)
% Step 3: Conduct Earth Gravity Assist (EGA) to lower orbit into a one-year-period
% Step 4: Perform one revolution around the sun, then burn to match Earth's velocity & rendezvous


%% HANDLE INPUTS
if nargin < 3
    % start_date not given - default is 1/1/2024
    start_date = datetime(2024,1,1);
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

% % Starting satellite orbit parameters (around Earth)
% switch parking_orbit
%     case "GTO"
%         parking.rp = earth_rad + 400; % Radius of perigee for parking orbit (km)
%         parking.ra = 42164; % Radius of apogee for parking orbit (km)
%     case "LEO"
%         parking.rp = earth_rad + 200; % Radius of perigee for parking orbit (km)
%         parking.ra = parking.rp; % Radius of apogee for parking orbit (km)
% end
% parking.sma = (parking.rp + parking.ra)/2; % semi-major axis (km)
% parking.ecc = (parking.ra - parking.rp)/(parking.ra + parking.rp); % eccentricity
% parking.eps = -earth_mu/(2*parking.sma); % parking orbital energy (km^2/s^2)


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


%% SUMMARY



%% Organize trajectory values (delta V's, important times in trajectory, orbital elements, etc. 
% for use in STK simulation









end




