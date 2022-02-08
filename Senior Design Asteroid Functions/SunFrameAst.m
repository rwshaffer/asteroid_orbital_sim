function [rT1, rT2, aT, ET, VT1_sun, VT2_sun] = SunFrameAst(a_E, a_A, a_trans, mu_S)
rT1 = a_E;
rT2 = a_A;
aT = a_trans;
ET = -mu_S/2/aT; %remember the body we are orbiting

VT1_sun = sqrt(2*(ET + mu_S/rT1)); %velocity of the transfer that the Sun sees (Earth)
VT2_sun = sqrt(2*(ET + mu_S/rT2)); %" " (Venus)
end