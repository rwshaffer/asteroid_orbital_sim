function [dV] = C3_to_dV(C3, init_orbit, r0, v0)

% Function to determine the impulsive delta-V needed to obtain a given C3 (characteristic energy) from a specified
% starting point (specified by either init_orbit or r0, v0)

% NOTE: Impulsive only right now

% init_orbit: Can be LEO, GEO, or Lunar. Do not specify r0 or v0 (only give two inputs)
% r0, v0 can be used to specify any arbitrary starting state. Input init_orbit as ""

% INPUTS:
% C3    -       Characteristic energy (km^2/s^2)
% init_orbit -  "LEO", "GEO", "Lunar", or ""
% r0    -       Initial radius (km)
% v0    -       Initial velocity (km/s)

% OUTPUTS:
% dV    -       Delta-V required (km/s)

muEarth = 3.986e5;
rEarth = 6378;

if nargin == 2

    switch init_orbit
        case "LEO"
            r0 = rEarth + 400;      
        case "GEO"
            r0 = rEarth + 35785;
        case "Lunar"
            r0 = 3.48e8;
    end
    v0 = sqrt(muEarth/r0);

elseif nargin ~= 4 || init_orbit ~= ""
    disp('Input Error')
    return
end

epsilon = C3/2; % Orbital energy for desired trajectory
v1 = sqrt(2*(epsilon + muEarth/r0));

dV = v1 - v0;


end