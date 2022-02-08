function [r3_A, v_circ_A, VT2_ast, r_SOI_ast, E_trans_ast, v3T_A, dV_retro] = AsteroidFrame(a_A, mu_A, VT2_sun, v_A, mu_S, vAtm, rA)
r3_A = vAtm + rA; %radius of the final parking orbit, L7 slide 22
v_circ_A = sqrt(mu_A/r3_A); %spacecraft circular velocity calculation
VT2_ast = VT2_sun - v_A; %getting the velocity from Sun's POV to Venus's; velocity at Venus SOI
r_SOI_ast = a_A*(mu_A/mu_S)^(2/5);
E_trans_ast = VT2_ast^2/2 - mu_A/r_SOI_ast; %km^2/s^2 

v3T_A = sqrt(2*((mu_A/r3_A) + E_trans_ast))
dV_retro = abs(v3T_A - v_circ_A); %instantaneous change in velocity, from transfer to circular after burn
% recall: deltaV = velocityNeededForTransfer - circularVelocitySCisGoing
fprintf('dV retro from Venus''s SOI is: %4.2f km/s \n', dV_retro)
end