%function [] = escape_to_earth_plane_MCS(root,sat)


% Note for EGA code: one of the inputs to this function will be the Earth escape date, 
% but we need to model starting from the EGA date, so just add one year to 
% the escape date as an approximation

ASTG = sat.Propagator;

ASTG.Options.DrawTrajectoryIn3D = 0; % Turn off drawing trajectory while calculating

MCS = ASTG.MainSequence;
MCS.RemoveAll;

%% ----------------------------------------------------------------------------------------------- %%
%% TARGET ASTEROID PLANE AT EARTH ESCAPE
% Insert and configure a target sequence that will vary the Earth escape elements
% such that once the S/C has left the Earth's SOI, it will have the same inclination
% and RAAN as the asteroid; therefore it is in the same orbital plane.

%% TARGET SEQUENCE: Targeting the Asteroid's Plane
ts_ast_plane = MCS.Insert('eVASegmentTypeTargetSequence','Target Asteroid Plane','-');


%% SET INITIAL STATE - Target Vector Outgoing Asymptote
% Insert inside Target Asteroid Plane.
initstate = ts_ast_plane.Segments.Insert('eVASegmentTypeInitialState','Escape Periapsis','-');
% Set initial state - Outgoing Asymptote elements.
initstate.SetElementType('eVAElementTypeTargetVectorOutgoingAsymptote');
initstate.OrbitEpoch = "7 Jul 2025 07:00:00.000";
initstate.Element.RadiusofPeriapsis = 6778.14;
initstate.Element.C3Energy = 5.01067;
initstate.Element.RAOutgoing = 3.41985;
initstate.Element.DeclinationOutgoing = -61.8156;
initstate.Element.VelocityAzimuthPeriapsis = 0;
initstate.Element.TrueAnomaly = 0;
% Set initial S/C parameters
initstate.SpacecraftParameters.DryMass = 224;
initstate.FuelTank.FuelMass = 170;
% Set control variables for target sequence
% Note: these sucked to find. In STK help search "eVAControlInitState" then find an "Enumeration" page
initstate.EnableControlParameter('eVAControlInitStateEpoch')
initstate.EnableControlParameter('eVAControlInitStateTargetVecOutC3')
initstate.EnableControlParameter('eVAControlInitStateTargetVecOutAsympRA')
initstate.EnableControlParameter('eVAControlInitStateTargetVecOutAsympDec')
initstate.EnableControlParameter('eVAControlInitStateTargetVecOutVelAzAtPeriapsis')


%% MANEUVER - Finite burn to edge of Earth's SOI
% Insert inside Target Asteroid Plane.
burn_to_EGA_exit = ts_ast_plane.Segments.Insert('eVASegmentTypeManeuver','Burn to EGA Exit','-');
% Change to Finite and set correct engine
burn_to_EGA_exit.SetManeuverType('eVAManeuverTypeFinite');
burn_to_EGA_exit.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
% Configure propagator and stopping conditions
burn_propagator = burn_to_EGA_exit.Maneuver.Propagator;
burn_propagator.MaxPropagationTime = 86400*100; % DONT FORGET THIS when using stopping conditions besides duration, epoch
burn_propagator.StoppingConditions.Add('R_Magnitude');
burn_propagator.StoppingConditions.Remove('Duration');
burn_propagator.StoppingConditions.Item('R_Magnitude').Properties.Trip = 1e6;


%% PROPAGATE - Brief propagation in Heliocentric orbit
% Insert inside Target Asteroid Plane.
exit_earth_SOI = ts_ast_plane.Segments.Insert('eVASegmentTypePropagate','Exit Earth SOI','-');
% Configure propagator and stopping condition
exit_earth_SOI.PropagatorName = 'Heliocentric';
exit_earth_SOI.StoppingConditions.Item('Duration').Properties.Trip = 60; % Very short propagate segment
% Configure results for target sequence - Inclination and RAAN in heliocentric orbit
exit_earth_SOI.Results.Add('Keplerian Elems/Inclination');
exit_earth_SOI.Results.Item('Inclination').CoordSystemName = 'CentralBody/Sun MeanEclpJ2000';
exit_earth_SOI.Results.Add('Keplerian Elems/RAAN');
exit_earth_SOI.Results.Item('RAAN').CoordSystemName = 'CentralBody/Sun MeanEclpJ2000';


%% Configure TS - Target Asteroid Plane
% Get handle to differential corrector used in target sequence
dc = ts_ast_plane.Profiles.Item('Differential Corrector');
% Configure Epoch control parameter - change its max step & perturbation
init_epoch_control = dc.ControlParameters.GetControlByPaths('Escape Periapsis','InitialState.Epoch');
init_epoch_control.Enable = true;
init_epoch_control.MaxStep = 3000;
init_epoch_control.Perturbation = 100;
% Save long string for control parameter names
init_control = 'InitialState.Target_Vector_Outgoing_Asymptote';
% Configure C3 control parameter - change its max step & perturbation
init_c3_control = dc.ControlParameters.GetControlByPaths('Escape Periapsis',sprintf('%s.C3',init_control));
init_c3_control.Enable = true;
init_c3_control.MaxStep = 0.5;
init_c3_control.Perturbation = 0.01;
% Enable other control parameters - no need to change max steps/perturbations
dc.ControlParameters.GetControlByPaths('Escape Periapsis',sprintf('%s.AsymDec',init_control)).Enable = true;
dc.ControlParameters.GetControlByPaths('Escape Periapsis',sprintf('%s.AsymRA',init_control)).Enable = true;
dc.ControlParameters.GetControlByPaths('Escape Periapsis',sprintf('%s.AzVp',init_control)).Enable = true;
% Enable results and configure so that desired values = ML values
inc_result = dc.Results.GetResultByPaths('Exit Earth SOI','Inclination');
inc_result.Enable = true;
inc_result.DesiredValue = ML.inc;
raan_result = dc.Results.GetResultByPaths('Exit Earth SOI','RAAN');
raan_result.Enable = true;
raan_result.DesiredValue = ML.raan;
% Set final DC and targeter properties and run modes
dc.MaxIterations = 50;
dc.EnableDisplayStatus = true;
dc.Mode = 'eVAProfileModeIterate';
ts_ast_plane.Action = 'eVATargetSeqActionRunActiveProfiles';


%% ----------------------------------------------------------------------------------------------- %%
%% TARGET ASTEROID ORBIT
% Insert and configure a target sequence that will vary the S/C thrust for roughly one orbit
% around the Sun such that at the S/C descending node (one orbit after Earth escape), 
% it will be very close to the asteroid target.
%% TARGET SEQUENCE: Targeting the Asteroid's Orbit
ts_ast_orbit = MCS.Insert('eVASegmentTypeTargetSequence','Target Asteroid Orbit','-');


%% MANEUVER - Finite burn until S/C apoapsis is near asteroid's orbit
% Insert inside Target Asteroid Orbit.
burn_raise_apo = ts_ast_orbit.Segments.Insert('eVASegmentTypeManeuver','Burn Raise Apoapsis','-');
% Change to Finite and set attitude thrust vector
burn_raise_apo.SetManeuverType('eVAManeuverTypeFinite');
burn_raise_apo.Maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
burn_raise_apo.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite VNC(Sun)';
%   Note: it appears that this defaults to spherical which is what we want. Unsure though
%   Note: need to find a way to set thrust vector values!
% Change to correct engine and configure thrust efficiency
burn_raise_apo.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_raise_apo.Maneuver.ThrustEfficiency = 0.8;
burn_raise_apo.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and stopping conditions - Use semimajor axis as a stopping condition
burn_propagator = burn_raise_apo.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Heliocentric';
burn_propagator.MaxPropagationTime = 86400*300;
burn_propagator.StoppingConditions.Add('UserSelect');
burn_propagator.StoppingConditions.Remove('Duration');
sma_stop = burn_propagator.StoppingConditions.Item('UserSelect');
sma_stop.Properties.UserCalcObjectName = 'Semimajor Axis';
sma_stop.Properties.UserCalcObject.CentralBodyName = 'Sun';
sma_stop.Properties.Trip = 1.78742e+08;
% Set control parameters for target sequence
burn_raise_apo.EnableControlParameter('eVAControlManeuverFiniteSphericalAz');
burn_raise_apo.EnableControlParameter('eVAControlManeuverFiniteSphericalElev');
burn_raise_apo.EnableControlParameter('eVAControlManeuverFiniteThrustEfficiency');
sma_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');


%% PROPAGATE - Coast to S/C Apoapsis
% Insert inside Target Asteroid Plane.
coast_to_apo = ts_ast_orbit.Segments.Insert('eVASegmentTypePropagate','Coast to Apoapsis','-');
% Configure propagator and stopping condition
coast_to_apo.PropagatorName = 'Heliocentric';
coast_to_apo.StoppingConditions.Add('Apoapsis');
coast_to_apo.StoppingConditions.Remove('Duration'); 
coast_to_apo.StoppingConditions.Item('Apoapsis').Properties.CentralBodyName = 'Sun';


%% MANEUVER - Finite burn to raise periapsis until reaching descending node
% Insert inside Target Asteroid Orbit.
burn_raise_peri = ts_ast_orbit.Segments.Insert('eVASegmentTypeManeuver','Burn Raise Periapsis','-');
% Change to Finite and set attitude thrust vector
burn_raise_peri.SetManeuverType('eVAManeuverTypeFinite');
burn_raise_peri.Maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
burn_raise_peri.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite VNC(Sun)';
%   Note: it appears that this defaults to spherical which is what we want. Unsure though
%   Note: need to find a way to set thrust vector values!
% Change to correct engine and configure thrust efficiency
burn_raise_peri.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_raise_peri.Maneuver.ThrustEfficiency = 0.1;
burn_raise_peri.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and stopping conditions - Use semimajor axis as a stopping condition
burn_propagator = burn_raise_peri.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Heliocentric';
burn_propagator.MaxPropagationTime = 86400*300;
burn_propagator.StoppingConditions.Add('DescendingNode');
burn_propagator.StoppingConditions.Remove('Duration');
burn_propagator.StoppingConditions.Item('DescendingNode').Properties.CoordSystem ...
    = 'CentralBody/Sun MeanEclpJ2000';
% Set control parameters for target sequence
burn_raise_peri.EnableControlParameter('eVAControlManeuverFiniteSphericalAz');
burn_raise_peri.EnableControlParameter('eVAControlManeuverFiniteSphericalElev');
burn_raise_peri.EnableControlParameter('eVAControlManeuverFiniteThrustEfficiency');
% Configure results for target sequence - Multibody and spherical elems for 1989 ML rendezvous
burn_raise_peri.Results.Add('MultiBody/Delta Declination');
burn_raise_peri.Results.Item('Delta_Declination').CentralBodyName = '1989_ML';
burn_raise_peri.Results.Add('MultiBody/Delta Right Asc');
burn_raise_peri.Results.Item('Delta_Right_Asc').CentralBodyName = '1989_ML';
burn_raise_peri.Results.Add('Spherical Elems/R Mag');
burn_raise_peri.Results.Item('R_Mag').ReferencePointName = 'CentralBody/1989_ML Center';
burn_raise_peri.Results.Add('Spherical Elems/V Mag');
burn_raise_peri.Results.Item('V_Mag').CoordSystemName = 'CentralBody/1989_ML J2000';


%% Configure TS - Target Asteroid Orbit
% Get handle to differential corrector used in target sequence
dc = ts_ast_orbit.Profiles.Item('Differential Corrector');
% Configure Thrust Efficiency control parameters - change their max steps & perturbations
burn1eff_control = dc.ControlParameters.GetControlByPaths('Burn_Raise_Apoapsis','FiniteMnvr.Thrusting.ThrustEfficiency');
burn1eff_control.Enable = true;
burn1eff_control.MaxStep = 0.005;
burn1eff_control.Perturbation = 0.0001;
burn2eff_control = dc.ControlParameters.GetControlByPaths('Burn_Raise_Periapsis','FiniteMnvr.Thrusting.ThrustEfficiency');
burn2eff_control.Enable = true;
burn2eff_control.MaxStep = 0.005;
burn2eff_control.Perturbation = 0.0001;
% Configure Semi-Major axis control parameter - change max step and perturbation
sma_control = dc.ControlParameters.GetControlByPaths('Burn_Raise_Apoapsis','FiniteMnvr.StoppingConditions.UserSelect.TripValue');
sma_control.Enable = true;
sma_control.MaxStep = 100000;
sma_control.Perturbation = 100;
% Enable thrust vector controls (no need to configure)
pointing_string = 'FiniteMnvr.Pointing.Spherical';
dc.ControlParameters.GetControlByPaths('Burn_Raise_Apoapsis',sprintf('%s.Azimuth',pointing_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_Raise_Apoapsis',sprintf('%s.Elevation',pointing_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_Raise_Periapsis',sprintf('%s.Azimuth',pointing_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_Raise_Periapsis',sprintf('%s.Elevation',pointing_string)).Enable = true;
% Enable results and configure tolerances
delta_dec_result = dc.Results.GetResultByPaths('Burn_Raise_Periapsis','Delta_Declination');
delta_dec_result.Enable = true;
delta_dec_result.Tolerance = 0.05;
delta_ra_result = dc.Results.GetResultByPaths('Burn_Raise_Periapsis','Delta_Right_Asc');
delta_ra_result.Enable = true;
delta_ra_result.Tolerance = 0.05;
rmag_result = dc.Results.GetResultByPaths('Burn_Raise_Periapsis','R_Mag');
rmag_result.Enable = true;
rmag_result.Tolerance = 1.5e5;
vmag_result = dc.Results.GetResultByPaths('Burn_Raise_Periapsis','V_Mag');
vmag_result.Enable = true;
vmag_result.Tolerance = 0.45;
% Set final DC and targeter properties and run modes
dc.MaxIterations = 50;
dc.EnableDisplayStatus = true;
dc.Mode = 'eVAProfileModeIterate';
ts_ast_orbit.Action = 'eVATargetSeqActionRunActiveProfiles';





ASTG.RunMCS

ASTG.RunMCS

%end