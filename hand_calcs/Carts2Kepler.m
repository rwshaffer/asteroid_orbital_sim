function [orbit,e_vec,h_vec,eps] = Carts2Kepler(R_vec,V_vec,central_body)
% Converts Cartesian position and velocity vectors to Keplerian elements
% in an inertial frame
%
% INPUTS
% R_vec: position vector
% V_vec: velocity vector
% central_body: 'earth' or 'sun'
%
% OUTPUTS
% orbit: struct containing Keplerian elements of orbit
%   orbit.epoch - epoch at which orbit is based (JD)
%   orbit.sma - semimajor axis (km)
%   orbit.ecc - eccentricity
%   orbit.inc - inclination (deg)
%   orbit.raan - right ascension of the ascending node (deg)
%   orbit.aop - argument of periapsis (deg)
%   orbit.nu - true anomaly (deg)
% e_vec: eccentricity vector (pointed to periapsis)
% h_vec: angular momentum vector (normal to trajectory)
% eps: orbital specific energy (km^2/s^2)

%% Helpful values and constants
r = norm(R_vec); % Magnitude of position
v = norm(V_vec); % Magnitude of velocity
i_hat = [1;0;0]; % Positive X-axis in ref. frame
k_hat = [0;0;1]; % Positive Z-axis in ref. frame
mu_earth = 3.986 * 10^5; % Earth Gravitational Constant
mu_sun =  1.3271233e11; % Sun Gravitational Constant

% Central body gravitational constant
switch central_body
    case 'sun'
        mu = mu_sun;
    case 'earth'
        mu = mu_earth;
end

%% Semi-major axis and eccentricity
h_vec = cross(R_vec,V_vec); % Angular momentum vector (normal to orbit)
h = norm(h_vec); % Angular momentum
eps = v^2/2 - mu/r; % Orbital energy (km^2/s^2)
sma = -mu/(2*eps); % Semi-major axis (km)
e_vec = 1/mu * ((v^2 - mu/r)*R_vec - (dot(R_vec,V_vec))*V_vec); % Eccentricity vector (pointed to periapsis)
ecc = norm(e_vec); % Eccentricity
%% True anomaly
nu = acos(dot(R_vec,e_vec)/(r*ecc)); % True anomaly (rad)
% Quadrant check
if dot(R_vec,V_vec) < 0
    nu = 2*pi - nu;
end

%% Inclination
inc = acos(dot(k_hat,h_vec)/h); % Inclination (rad)
%% Right Ascension of the Ascending Node
n_vec = cross(k_hat,h_vec); % intersection of orbital plane and reference XY plane
raan = acos(dot(i_hat,n_vec)/norm(n_vec)); % Right Ascension (rad)
% Quadrant check
if n_vec(2) < 0
    raan = 2*pi - raan;
end
    
%% Argument of Periapsis
aop = acos(dot(e_vec,n_vec)/(ecc*norm(n_vec))); % Argument of periapsis (rad)
% Quadrant check
if e_vec(3) < 0
    aop = 2*pi - aop;
end

%% Pack elements into output struct
orbit.sma = sma;
orbit.ecc = ecc;
orbit.inc = inc;
orbit.raan = raan;
orbit.aop = aop;
orbit.nu = nu;


end