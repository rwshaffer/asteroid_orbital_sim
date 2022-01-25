function [TOF] = TimeOfFlight(orbit,mu,final_nu,K)
% Given an orbit, an initial and a final true anomaly, 
%   what is the time of flight?
%
% INPUTS:
% orbit: struct containing Keplerian elements of orbit
%   orbit.epoch - epoch at which orbit is based
%   orbit.sma - semimajor axis (km)
%   orbit.ecc - eccentricity
%   orbit.inc - inclination (rad)
%   orbit.raan - right ascension of the ascending node (rad)
%   orbit.aop - argument of periapsis (rad)
%   orbit.nu - initial true anomaly (rad)
% mu: gravitational parameter (km^3/s^2)
% final_nu: final true anomaly (rad)
% K: number of times that the object passes periapsis. Default 0.
%
% OUTPUTS:
% TOF: time of flight (seconds)

if nargin < 4
    K = 0;
end

[~,~,M_0] = convert_anomalies(orbit.nu,orbit.ecc,"mean"); % Initial mean anomaly (rad)
[~,~,M_f] = convert_anomalies(final_nu,orbit.ecc,"mean"); % Final mean anomaly (rad)
TOF = sqrt(orbit.sma^3/mu) * (2*pi*K + M_f - M_0); % Time of flight (sec)


end