clear; close all; clc


load Earth_Orbit_circ.mat
load ML_Orbit.mat
target = ML;

% Pick a date for the EGA
% Bunch of Lambert Problems to get from Earth to 1989 ML
% See what delta-V is needed to match the Earth rel. velocity magnitude
% Find a trajectory with a period of 365 days with that rel. velocity

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



%% EGA 
ega_date = datetime(2025,7,1); % Test
ega_jd = juliandate(ega_date);

post_ega_TOFs = linspace(50,730,100);

%% Propagate Earth to position at EGA
earth.ega_tof = (ega_jd - earth.epoch)*86400; % Time (sec) between Earth epoch and EGA
earth.nu_ega = KeplerProblem(earth,sun_mu,earth.ega_tof); % Target true anomaly at launch time (rad)
earth.nu = earth.nu_ega; % Update true anomaly to position at launch
earth.epoch = ega_jd; % Update target epoch to launch window

%% Iterate through Lambert Problem times of flight

% Save target's epoch conditions
target.init_epoch = target.epoch;
target.init_nu = target.nu;

% Initialize vectors for storing outputs
rel_earth_velocities = zeros(3,length(post_ega_TOFS));
rdvs_delta_vs = zeros(size(post_ega_TOFs));

for i = 1:length(post_ega_TOFs)
    % Propagate target to position at end of EGA
    traj_end_time = ega_jd + post_ega_TOFs(i); % JD at which S/C intercepts target
    target.prop_time = (traj_end_time - target.epoch)*86400; % Time (sec) between Earth epoch and EGA
    target.nu_final = KeplerProblem(target,sun_mu,target.prop_time); % Target true anomaly at launch time (rad)
    target.nu = target.nu_final; % Update true anomaly to position at launch
    target.epoch = traj_end_time; % Update target epoch to launch window

    % Convert orbital elements into pos/vel for Earth and target
    [earth.ega_pos,earth.ega_vel] = Kepler2Carts(earth,'sun');
    [target.final_pos,target.final_vel] = Kepler2Carts(target,'sun');
    
    % Lambert Problem to obtain velocities on transfer trajectory
    [V1,V2]=GaussProblemtextbook(earth.ega_pos,target.final_pos,post_ega_TOFs(i));

    % Delta V at asteroid rendezvous
    satellite.rdvs_delta_v = norm(V2 - target.final_vel);
    
    % Velocity relative to earth after EGA
    satellite.vel_rel_earth = V1 - earth.ega_vel;
    
    % Save relevant Lambert problem quantities in matrices
    rel_earth_velocities(:,i) = satellite.vel_rel_earth;
    rdvs_delta_vs(i) = satellite.rdvs_delta_v;
end

%% Find trajectories from Earth to Earth with each relative Earth velocity

% Velocity at a set radius as a function of orbital period
TP = 365 * 86400; % Orbital period of Earth to EGA transfer orbit (sec)
v = sqrt((-4*pi^2*sun_mu^2/TP^2)^(-1/3)+2*sun_mu/earth.ega_pos);

for i = 1:length(rel_earth_velocities)
    
    
    
    
end

