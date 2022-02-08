function [V_at_orbit, eFPA, V_at_EI] = Ast_EDL_InitialConditions(mu_A, E_trans_ast, r3_A, rEI)
V_at_orbit = sqrt(2*((mu_A/rEI) + E_trans_ast)); %velocity of the spacecraft as it enters the EI rel to the planet
% coming off hyperbolic trajectory ^^
% 10 km/s of delta v is a lot; will create a plasma where shielding is
% necessary. heat shield will cost extra mass, but we should not design to
% slow down anyway as fuel will cost more weight and money...
a_hyp = -mu_A/(2*E_trans_ast); %rearranged energy equation to get hyperbolic semimajor axis
e_hyp = 1-(r3_A/a_hyp);
p = a_hyp*(1-e_hyp^2);

h_hyp = sqrt(mu_A*p);

V_at_EI = sqrt(2*(E_trans_ast + (mu_A/rEI)));

eFPA = -acosd(h_hyp/(V_at_EI*(rEI))); %for some reason what is inside this is 1, causing eFPA at any alt to be 0 deg
fprintf('The velocity of the spacecraft at the entry interface is: %4.2f km/s \n', V_at_orbit)
fprintf('Based on these conditions, the FPA at the entry interface is: %4.2f degrees \n', eFPA)
end