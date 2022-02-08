%% E2AConstantsFunc
%% Sun
mu_S = 132712E6; %km^3/s^2

%% Earth
a_E = 149.6E6; %km
v_E = sqrt(mu_S/a_E); %km/s
mu_E = .3986E6; %km^3/s^2
w_E = sqrt(mu_S/a_E^3); %rad/s

%% Asteroid...CHANGE based on researched values
%These are Venus parameters...
%WHEN POSSIBLE, DO NOT HARDCODE IN VALUES
rA = 6051.8; %km, found on the internet
a_A = 108.2E6; %km
v_A = sqrt(mu_S/a_A); %km/s
mu_A = .32456E6; %km^3/s^2
w_A = sqrt(mu_S/a_A^3); %rad/s