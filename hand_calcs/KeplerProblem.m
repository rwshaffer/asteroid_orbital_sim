function [nu_final] = KeplerProblem(orbit,mu,TOF)
% Solve the Kepler Problem: given orbital elements (including true anomaly)
% at an initial time, and a time of flight, what is the final true anomaly?
% Used to propagate an orbit through some time interval.

% INPUTS
% orbit: struct containing Keplerian elements of orbit
%   orbit.epoch - epoch at which orbit is based
%   orbit.sma - semimajor axis (km)
%   orbit.ecc - eccentricity
%   orbit.inc - inclination (deg)
%   orbit.raan - right ascension of the ascending node (deg)
%   orbit.aop - argument of periapsis (deg)
%   orbit.nu - true anomaly (deg)
% mu: Gravitational parameter of central body (km^3/s^2)
% TOF: time of flight (sec)
%
% OUTPUTS
% nu_final: final true anomaly (deg)

%% Find Mean Anomaly (M)
n_bar = sqrt(mu/orbit.sma^3); % Mean motion (rad/sec)
delta_M = n_bar * TOF; % Change in Mean Anomaly (rad)

%% Find new True Anomaly (nu)
[delta_nu,~,~] = convert_anomalies(delta_M,orbit.ecc,"mean"); % Change in true anomaly (rad)
delta_nu = delta_nu * 180/pi; % Change in true anomaly (deg)
nu_final = mod(orbit.nu + delta_nu,360); % Final true anomaly (deg)
end