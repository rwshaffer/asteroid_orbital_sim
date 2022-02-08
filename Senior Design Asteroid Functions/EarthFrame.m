function [v_circ_E, VT1_earth, r_SOI_E, E_trans_Earth, v0T_E, dV_boost] = EarthFrame(mu_E, VT1_sun, v_E, a_E, mu_S, alt)
r0_E = alt + 6378; %estimation of satellite's altitude (in LEO) 
v_circ_E = sqrt(mu_E/r0_E);
VT1_earth = VT1_sun - v_E; %getting the velocity from Sun's POV to Earth's
r_SOI_E = a_E*(mu_E/mu_S)^(2/5);
E_trans_Earth = VT1_earth^2/2 - mu_E/r_SOI_E ; %km^2/s^2 

v0T_E = sqrt(2*((mu_E/r0_E) + E_trans_Earth));
dV_boost = abs(v0T_E - v_circ_E); %instantaneous change in velocity, from circular to transfer after boost
fprintf('dV boost from the Earth''s SOI is: %4.2f km/s \n', dV_boost)
