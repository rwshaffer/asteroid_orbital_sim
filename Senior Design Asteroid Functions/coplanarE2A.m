function [TOF, alpha_lead_d, phi_impulse, T_wait_year, a_trans] = coplanarE2A(a_E, a_A, mu_S, w_A, w_E)
a_trans = (a_E + a_A)/2;
TOF = pi*sqrt(a_trans^3/mu_S);
fprintf('The Time of Flight of this problem is: %4.2f days \n', TOF/3600/24)
%---------------------
w_target = w_A;
alpha_lead = w_target*TOF;
alpha_lead_d = (180/pi)*(w_target*TOF);
fprintf('The lead angle of the target is: %4.2f radians, or %8.3f degrees \n', alpha_lead, alpha_lead_d)
%---------------------
phi_impulse = pi - alpha_lead;
fprintf('phi_impulse is: %4.2f radians \n', phi_impulse) %how do i change decimal precision for fprintf? 
%---------------------
phi_init = 0; %rad, beginning angle
w_interceptor = w_E;
T_wait = (phi_impulse - phi_init + 2*pi)/(w_target - w_interceptor);
T_wait_year = T_wait/(3600*24*365); 
fprintf('wait time is: %4.2f seconds, which is %8.3f years \n', T_wait, T_wait_year)
end
