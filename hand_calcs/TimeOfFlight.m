function [TOF] = TimeOfFlight(orbit,mu,final_nu,K)
% Given an orbit, an initial and a final true anomaly, 
%   what is the time of flight?
%
% INPUTS:
% orbit: struct containing Keplerian elements of orbit
%   orbit.epoch - epoch at which orbit is based
%   orbit.sma - semimajor axis (km)
%   orbit.ecc - eccentricity
%   orbit.inc - inclination (deg)
%   orbit.raan - right ascension of the ascending node (deg)
%   orbit.aop - argument of periapsis (deg)
%   orbit.nu - initial true anomaly (deg)
% mu: gravitational parameter (km^3/s^2)
% final_nu: final true anomaly (deg)
% K: number of times that the object passes periapsis. Default 0.
%
% OUTPUTS:
% TOF: time of flight (seconds)

if nargin < 4
    K = 0;
end

nu_0 = orbit.nu * pi/180; % Initial true anomaly (rad)
[~,~,M_0] = convert_anomalies(nu_0,orbit.ecc,"mean"); % Initial mean anomaly (rad)
nu_f = final_nu * pi/180; % Final true anomaly (rad)
[~,~,M_f] = convert_anomalies(nu_f,orbit.ecc,"mean"); % Final mean anomaly (rad)
TOF = sqrt(orbit.sma^3/mu) * (2*pi*K + M_f - M_0); % Time of flight (sec)


end