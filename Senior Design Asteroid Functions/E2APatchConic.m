%% The patch-conic problem, from Earth to Asteroid
%Genevieve Gemond, Virginia Tech Aerospace Engineering Class of 2022
clear; home
%% Sun, Earth, Venus constants 
E2AConstantsFunc
alt0 = 600; %km, initial parking orbit altitude about Earth. Changes to this will have minor changes in final results. 
aAtm = 200; %km, CHANGE thickness of Venus atmosphere; the desired entry interface
%% Part 1: Coplanar rendezvous *Around the Sun*
%plan trajectory
[TOF, alpha_lead_d, phi_impulse, T_wait_year, a_trans] = coplanarE2A(a_E, a_A, mu_S, w_A, w_E); 
%coplanar to achieve hohmann, optimized at 180 deg apart
%coplanar tells angular separation
%----
%execute
[rT1, rT2, aT, ET, VT1_sun, VT2_sun] = SunFrameAst(a_E, a_A, a_trans, mu_S);
%% Part 2a: 
%#1. Earth
%2a is doing everything relative to planet I am departing
%most likely no need to change any of this, but check if you would like 
[v_circ_E, VT1_earth, r_SOI_E, E_trans_Earth, v0T_E, dV_boost] = EarthFrame(mu_E, VT1_sun, v_E, a_E, mu_S, alt0);
%% Part 2b: Zooming in on Individual Frames
%doing everything relative to planet I am arriving at
%zoom in on ASTEROID frame...CHANGE 
% #2. Venus
[r3_A, v_circ_A, VT2_ast, r_SOI_ast, E_trans_ast, v3T_A, dV_retro] = AsteroidFrame(a_A, mu_A, VT2_sun, v_A, mu_S, aAtm, rA)

%% Part 2 - getting initial contidions for EDL simulation 
%Velocity of Spacecraft @ Venus's SOI
%If we want EI conditions, CHANGE 

rEI = r3_A + aAtm;
[V_at_orbit, eFPA, V_at_EI] = Ast_EDL_InitialConditions(mu_A, E_trans_ast, r3_A, rEI);