function [sat,diverged] = escape_to_earth_plane_MCS(ML,sat,launch_date)


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
initstate.OrbitEpoch = launch_date;
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
%burn_propagator.MaxPropagationTime = 86400*100; % DONT FORGET THIS when using stopping conditions besides duration, epoch
burn_propagator.StoppingConditions.Add('R_Magnitude');
burn_propagator.StoppingConditions.Remove('Duration');
burn_propagator.StoppingConditions.Item('R_Magnitude').Properties.Trip = 1e6;
burn_propagator.Results.Add('Maneuver/DeltaV');


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
burn_raise_apo.Maneuver.SetAttitudeControlType('eVAAttitudeControlTimeVarying');
burn_raise_apo.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite VNC(Sun)';
burn_raise_apo.Maneuver.AttitudeControl.Az0 = -4.484;
burn_raise_apo.Maneuver.AttitudeControl.El0 = -14.557;

%   Note: it appears that this defaults to spherical which is what we want. Unsure though
% Change to correct engine and configure thrust efficiency
burn_raise_apo.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_raise_apo.Maneuver.ThrustEfficiency = 0.776;
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
burn_raise_apo.Results.Add('Maneuver/DeltaV');
% Set control parameters for target sequence
burn_raise_apo.EnableControlParameter('eVAControlManeuverFiniteAz0');
burn_raise_apo.EnableControlParameter('eVAControlManeuverFiniteEl0');
burn_raise_apo.EnableControlParameter('eVAControlManeuverFiniteThrustEfficiency');
sma_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');


%% PROPAGATE - Coast to S/C Apoapsis
% Insert inside Target Asteroid Plane.
coast_to_apo = ts_ast_orbit.Segments.Insert('eVASegmentTypePropagate','Coast to Apoapsis','-');
% Configure propagator and stopping condition
coast_to_apo.PropagatorName = 'Heliocentric';
coast_to_apo.MaxPropagationTime = 86400*300;
coast_to_apo.StoppingConditions.Add('Apoapsis');
coast_to_apo.StoppingConditions.Remove('Duration'); 
coast_to_apo.StoppingConditions.Item('Apoapsis').Properties.CentralBodyName = 'Sun';


%% MANEUVER - Finite burn to raise periapsis until reaching descending node
% Insert inside Target Asteroid Orbit.
burn_raise_peri = ts_ast_orbit.Segments.Insert('eVASegmentTypeManeuver','Burn Raise Periapsis','-');
% Change to Finite and set attitude thrust vector
burn_raise_peri.SetManeuverType('eVAManeuverTypeFinite');
burn_raise_peri.Maneuver.SetAttitudeControlType('eVAAttitudeControlTimeVarying');
burn_raise_peri.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite VNC(Sun)';
burn_raise_peri.Maneuver.AttitudeControl.Az0 = -10.99;
burn_raise_peri.Maneuver.AttitudeControl.El0 = 0.901;
% Change to correct engine and configure thrust efficiency
burn_raise_peri.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_raise_peri.Maneuver.ThrustEfficiency = 0.108;
burn_raise_peri.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and stopping conditions - Use descending node as a stopping condition
burn_propagator = burn_raise_peri.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Heliocentric';
burn_propagator.MaxPropagationTime = 86400*300;
burn_propagator.StoppingConditions.Add('DescendingNode');
burn_propagator.StoppingConditions.Remove('Duration');
burn_propagator.StoppingConditions.Item('DescendingNode').Properties.CoordSystem ...
    = 'CentralBody/Sun MeanEclpJ2000';
% Set control parameters for target sequence
burn_raise_peri.EnableControlParameter('eVAControlManeuverFiniteAz0');
burn_raise_peri.EnableControlParameter('eVAControlManeuverFiniteEl0');
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
burn_raise_peri.Results.Add('Maneuver/DeltaV');

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
pointing_string = 'FiniteMnvr.Pointing.TimeVarying';
dc.ControlParameters.GetControlByPaths('Burn_Raise_Apoapsis',sprintf('%s.Az0',pointing_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_Raise_Apoapsis',sprintf('%s.El0',pointing_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_Raise_Periapsis',sprintf('%s.Az0',pointing_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_Raise_Periapsis',sprintf('%s.El0',pointing_string)).Enable = true;
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


%% ----------------------------------------------------------------------------------------------- %%
%% TARGET RENDEZVOUUS
% Insert and configure a target sequence that will vary the impulsive and finite burns
% when approaching the asteroid such that the orbital period roughly matches the asteroid, 
% and the relative distance and velocity are quite small (but still non-zero)
%% TARGET SEQUENCE: Targeting Rendezvous
ts_rdvs = MCS.Insert('eVASegmentTypeTargetSequence','Target Rendezvous','-');


%% MANEUVER - Impulsive burn
% Insert inside Target Asteroid Orbit.
impulsive_burn = ts_rdvs.Segments.Insert('eVASegmentTypeManeuver','Impulsive Burn','-');
% Set attitude thrust vector and magnitude
impulsive_burn.Maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
impulsive_burn.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite/AU_Digger_Sat VNC(1989_ML)';
impulsive_burn.Maneuver.AttitudeControl.CoordType = 'eVASphericalImpDeltaV';
impulsive_burn.Maneuver.AttitudeControl.Azimuth = 207.7;
impulsive_burn.Maneuver.AttitudeControl.Elevation = 26.8373;
impulsive_burn.Maneuver.AttitudeControl.Magnitude = 190;
impulsive_burn.Maneuver.AttitudeControl.AllowNegativeSphericalMagnitude = 1;
impulsive_burn.Results.Add('Maneuver/DeltaV');
% Configure control parameters for target sequence
impulsive_burn.EnableControlParameter('eVAControlManeuverImpulsiveSphericalAz');
impulsive_burn.EnableControlParameter('eVAControlManeuverImpulsiveSphericalElev');
impulsive_burn.EnableControlParameter('eVAControlManeuverImpulsiveSphericalMag');


%% MANEUVER - Finite burn to "rendezvous" with asteriod (closest we get in interplanetary traj.)
% Insert inside Target Rendezvous.
burn_close_app = ts_rdvs.Segments.Insert('eVASegmentTypeManeuver','Burn Close Approach','-');
% Change to Finite and set attitude thrust vector
burn_close_app.SetManeuverType('eVAManeuverTypeFinite');
burn_close_app.Maneuver.SetAttitudeControlType('eVAAttitudeControlTimeVarying');
burn_close_app.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite/AU_Digger_Sat VNC(1989_ML)';
burn_close_app.Maneuver.AttitudeControl.Az0 = 153.255;
burn_close_app.Maneuver.AttitudeControl.El0 = 6.83;
%burn_close_app.Maneuver.AttitudeControl.
% Change to correct engine and configure thrust efficiency
burn_close_app.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_close_app.Maneuver.ThrustEfficiency = 0.5;
burn_close_app.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and stopping conditions - Use Epoch as a stopping condition
burn_propagator = burn_close_app.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Heliocentric';
burn_propagator.MaxPropagationTime = 86400*300;
burn_propagator.StoppingConditions.Add('Epoch');
burn_propagator.StoppingConditions.Remove('Duration');
epoch_stop = burn_propagator.StoppingConditions.Item('Epoch');
epoch_stop.Properties.Trip = "17 Nov 2026 07:43:00.000";
% Set control parameters for target sequence
burn_close_app.EnableControlParameter('eVAControlManeuverFiniteAz0');
burn_close_app.EnableControlParameter('eVAControlManeuverFiniteEl0');
epoch_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');

% Configure results for target sequence - Multibody and spherical elems for 1989 ML rendezvous
burn_close_app.Results.Add('Keplerian Elems/Orbit_Period');
burn_close_app.Results.Item('Orbit_Period').CentralBodyName = 'Sun';
burn_close_app.Results.Add('Spherical Elems/R Mag');
burn_close_app.Results.Item('R_Mag').ReferencePointName = 'CentralBody/1989_ML Center';
burn_close_app.Results.Add('Spherical Elems/V Mag');
burn_close_app.Results.Item('V_Mag').CoordSystemName = 'CentralBody/1989_ML J2000';
burn_close_app.Results.Add('Maneuver/DeltaV');


%% Configure TS - Target Rendezvous
% Get handle to differential corrector used in target sequence
dc = ts_rdvs.Profiles.Item('Differential Corrector');
% Configure final Epoch control parameter - change max step
epoch_control = dc.ControlParameters.GetControlByPaths('Burn_Close_Approach','FiniteMnvr.StoppingConditions.Epoch.TripValue');
epoch_control.Enable = true;
epoch_control.MaxStep = 5000;
% Enable impulsive maneuver controls (no special configurations)
impulsive_string = 'ImpulsiveMnvr.Pointing.Spherical';
dc.ControlParameters.GetControlByPaths('Impulsive_Burn',sprintf('%s.Azimuth',impulsive_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Impulsive_Burn',sprintf('%s.Elevation',impulsive_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Impulsive_Burn',sprintf('%s.Magnitude',impulsive_string)).Enable = true;
% Enable finite thrust vector controls (no need to configure)
pointing_string = 'FiniteMnvr.Pointing.TimeVarying';
dc.ControlParameters.GetControlByPaths('Burn_Close_Approach',sprintf('%s.Az0',pointing_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_Close_Approach',sprintf('%s.El0',pointing_string)).Enable = true;
% Enable results and configure tolerances
period_result = dc.Results.GetResultByPaths('Burn_Close_Approach','Orbit_Period');
period_result.Enable = true;
period_result.DesiredValue = ML.period; % Note: in seconds
period_result.Tolerance = 0.75 * 86400; % 0.75 days, converted to seconds
rmag_result = dc.Results.GetResultByPaths('Burn_Close_Approach','R_Mag');
rmag_result.Enable = true;
rmag_result.Tolerance = 15400;
vmag_result = dc.Results.GetResultByPaths('Burn_Close_Approach','V_Mag');
vmag_result.Enable = true;
vmag_result.Tolerance = 0.1;
% Set final DC and targeter properties and run modes
dc.MaxIterations = 50;
dc.EnableDisplayStatus = true;
dc.Mode = 'eVAProfileModeIterate';


%% ----------------------------------------------------------------------------------------------- %%
%% TARGET ASTEROID PERIOD
% Insert and configure a target sequence that will vary the impulsive burn
% after closest arrival such that the orbital period exactly matches the asteroid
%% TARGET SEQUENCE: Targeting Period
ts_period = MCS.Insert('eVASegmentTypeTargetSequence','Target Orbital Period','-');


%% MANEUVER - Impulsive burn
% Insert inside Target Asteroid Orbit.
impulsive_burn_period = ts_period.Segments.Insert('eVASegmentTypeManeuver','Impulsive Burn Match Period','-');
% Set attitude thrust vector and magnitude
impulsive_burn_period.Maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
impulsive_burn_period.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite/AU_Digger_Sat VNC(1989_ML)';
impulsive_burn_period.Maneuver.AttitudeControl.CoordType = 'eVACartesianImpDeltaV';
impulsive_burn_period.Maneuver.AttitudeControl.X = -24.9;
impulsive_burn_period.Maneuver.AttitudeControl.Y = 2.51;
impulsive_burn_period.Maneuver.AttitudeControl.Z = 4.74;
% Configure control parameters for target sequence
impulsive_burn_period.EnableControlParameter('eVAControlManeuverImpulsiveCartesianX');
impulsive_burn_period.EnableControlParameter('eVAControlManeuverImpulsiveCartesianY');
impulsive_burn_period.EnableControlParameter('eVAControlManeuverImpulsiveCartesianZ');
% Configure results for target sequence - Period around sun and speed rel. to asteroid
impulsive_burn_period.Results.Add('Keplerian Elems/Orbit_Period');
impulsive_burn_period.Results.Item('Orbit_Period').CentralBodyName = 'Sun';
impulsive_burn_period.Results.Add('Spherical Elems/V Mag');
impulsive_burn_period.Results.Item('V_Mag').CoordSystemName = 'CentralBody/1989_ML J2000';
impulsive_burn_period.Results.Add('Maneuver/DeltaV');


%% Configure TS - Target Rendezvous
% Get handle to differential corrector used in target sequence
dc = ts_period.Profiles.Item('Differential Corrector');
% Enable impulsive maneuver controls and decrease perturbations
impulsive_string = 'ImpulsiveMnvr.Pointing.Cartesian';
x_control = dc.ControlParameters.GetControlByPaths('Impulsive_Burn_Match_Period',sprintf('%s.X',impulsive_string));
x_control.Enable = true;
x_control.Perturbation = 0.01;
y_control = dc.ControlParameters.GetControlByPaths('Impulsive_Burn_Match_Period',sprintf('%s.Y',impulsive_string));
y_control.Enable = true;
y_control.Perturbation = 0.01;
z_control = dc.ControlParameters.GetControlByPaths('Impulsive_Burn_Match_Period',sprintf('%s.Z',impulsive_string));
z_control.Enable = true;
z_control. Perturbation = 0.01;
% Enable results and configure tolerances
period_result = dc.Results.GetResultByPaths('Impulsive_Burn_Match_Period','Orbit_Period');
period_result.Enable = true;
period_result.DesiredValue = ML.period; % Note: in seconds
period_result.Tolerance = 10; % Match 1989 ML orbital period to within 10 seconds
vmag_result = dc.Results.GetResultByPaths('Impulsive_Burn_Match_Period','V_Mag');
vmag_result.Enable = true;
vmag_result.Tolerance = 0.005; % Note: in km/s
% Set final DC and targeter properties and run modes
dc.MaxIterations = 50;
dc.EnableDisplayStatus = true;
dc.Mode = 'eVAProfileModeIterate';


%% ----------------------------------------------------------------------------------------------- %%
%% TARGET ASTEROID PROX OPS
% Insert and configure a target sequence that will vary the impulsive burn
% after closest arrival such that the satellite velocity exactly matches the asteroid
%% TARGET SEQUENCE: Targeting Period
ts_proxops = MCS.Insert('eVASegmentTypeTargetSequence','Target Prox Ops','-');


%% MANEUVER - Impulsive burn
% Insert inside Target Asteroid Orbit.
impulsive_burn_match_speed = ts_proxops.Segments.Insert('eVASegmentTypeManeuver','Impulsive Burn Match Speed','-');
% Set attitude thrust vector and magnitude
impulsive_burn_match_speed.Maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
impulsive_burn_match_speed.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite/AU_Digger_Sat VNC(1989_ML)';
impulsive_burn_match_speed.Maneuver.AttitudeControl.CoordType = 'eVACartesianImpDeltaV';
impulsive_burn_match_speed.Maneuver.AttitudeControl.X = -7.035;
impulsive_burn_match_speed.Maneuver.AttitudeControl.Y = -0.05;
impulsive_burn_match_speed.Maneuver.AttitudeControl.Z = -0.05;
% Configure control parameters for target sequence
impulsive_burn_match_speed.EnableControlParameter('eVAControlManeuverImpulsiveCartesianX');
impulsive_burn_match_speed.EnableControlParameter('eVAControlManeuverImpulsiveCartesianY');
impulsive_burn_match_speed.EnableControlParameter('eVAControlManeuverImpulsiveCartesianZ');
% Configure results for target sequence - Period around sun and speed rel. to asteroid
impulsive_burn_match_speed.Results.Add('Spherical Elems/V Mag');
impulsive_burn_match_speed.Results.Item('V_Mag').CoordSystemName = 'CentralBody/1989_ML J2000';
impulsive_burn_match_speed.Results.Add('Maneuver/DeltaV');


%% Configure TS - Target Rendezvous
% Get handle to differential corrector used in target sequence
dc = ts_proxops.Profiles.Item('Differential Corrector');
% Enable impulsive maneuver controls (no special configurations)
impulsive_string = 'ImpulsiveMnvr.Pointing.Cartesian';
dc.ControlParameters.GetControlByPaths('Impulsive_Burn_Match_Speed',sprintf('%s.Z',impulsive_string)).Enable = true;
x_control = dc.ControlParameters.GetControlByPaths('Impulsive_Burn_Match_Speed',sprintf('%s.X',impulsive_string));
x_control.Enable = true;
x_control.Perturbation = 0.001;
y_control = dc.ControlParameters.GetControlByPaths('Impulsive_Burn_Match_Speed',sprintf('%s.Y',impulsive_string));
y_control.Enable = true;
y_control.Perturbation = 0.001;
z_control = dc.ControlParameters.GetControlByPaths('Impulsive_Burn_Match_Speed',sprintf('%s.Z',impulsive_string));
z_control.Enable = true;
z_control. Perturbation = 0.001;


% Enable results and configure tolerances
velocity_result = dc.Results.GetResultByPaths('Impulsive_Burn_Match_Speed','V_Mag');
velocity_result.Enable = true;
velocity_result.Tolerance = 3e-3;
% Set final DC and targeter properties and run modes
dc.MaxIterations = 50;
dc.EnableDisplayStatus = true;
dc.Mode = 'eVAProfileModeIterate';




%% Call function to systematically achieve convergence with all target sequences

[sat,diverged] = run_target_sequences(sat);




end