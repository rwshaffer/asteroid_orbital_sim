function [sat,diverged] = grav_assist_MCS(ML,sat,launch_date)

%% Run escape_to_earth_plane_MCS to get baseline trajectory

[sat,diverged] = escape_to_earth_plane_MCS(ML,sat,launch_date);

if diverged
    return
end
    
ASTG = sat.Propagator;

ASTG.Options.DrawTrajectoryIn3D = 0; % Turn off drawing trajectory while calculating
ASTG.Options.SmartRunMOde = 'eVASmartRunModeOnlyChanged'; % Run entire MCS


MCS = ASTG.MainSequence;

%% BACKWARD SEQUENCE
bs_ega = MCS.Insert('eVASegmentTypeBackwardSequence','Backwards from EGA','-');

%% TARGET SEQUENCE: Target Earth Escape from EGA
ts_escape = bs_ega.Segments.Insert('eVASegmentTypeTargetSequence','Target Earth Escape','-');

%% INITIAL STATE
backwards_init = ts_escape.Segments.Insert('eVASegmentTypeInitialState','Spawn at EGA','-');
% Find handle on forwards sequence initial state
initstate = MCS.Item(0).Segments.Item(0);
% Set initial state of backwards seq. to exactly match the forwards initial state
backwards_init.SetElementType('eVAElementTypeTargetVectorOutgoingAsymptote');
backwards_init.OrbitEpoch = launch_date;
backwards_init.Element.RadiusofPeriapsis = initstate.Element.RadiusofPeriapsis;
backwards_init.Element.C3Energy = initstate.Element.C3Energy;
backwards_init.Element.RAOutgoing = initstate.Element.RAOutgoing;
backwards_init.Element.DeclinationOutgoing = initstate.Element.DeclinationOutgoing;
backwards_init.Element.VelocityAzimuthPeriapsis = initstate.Element.VelocityAzimuthPeriapsis;
backwards_init.Element.TrueAnomaly = initstate.Element.TrueAnomaly;
backwards_init.SpacecraftParameters.DryMass = initstate.SpacecraftParameters.DryMass;
backwards_init.FuelTank.FuelMass = initstate.FuelTank.FuelMass;


%% IMPULSIVE MANEUVER AT EGA PERIAPSIS
dv_periapsis = ts_escape.Segments.Insert('eVASegmentTypeManeuver','Impulsive at Periapsis','-');
% Set magnitude (leave attitude control as "along thrust vector")
dv_periapsis.Maneuver.AttitudeControl.DeltaVMagnitude = .2286;
dv_periapsis.Results.Add('Maneuver/DeltaV');
% Configure control parameters for target sequence
dv_periapsis.EnableControlParameter('eVAControlManeuverImpulsiveSphericalMag');

%% FINITE MANEUVER TO EGA ENTRY
burn_ega_entry = ts_escape.Segments.Insert('eVASegmentTypeManeuver','Burn to EGA Entry','-');
% Change to Finite and leave attitude control as velocity vector
burn_ega_entry.SetManeuverType('eVAManeuverTypeFinite');
% Change to correct engine and configure thrust efficiency
burn_ega_entry.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_ega_entry.Maneuver.ThrustEfficiency = 1;
burn_ega_entry.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and stopping conditions - Use duration as a stopping condition
burn_propagator = burn_ega_entry.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Earth HPOP Default v10';
burn_propagator.StoppingConditions.Add('R_Magnitude');
burn_propagator.StoppingConditions.Remove('Duration');
burn_propagator.StoppingConditions.Item('R_Magnitude').Properties.Trip = 1e6;
burn_ega_entry.Results.Add('Maneuver/DeltaV');

%% FINITE MANEUVER TO APOAPSIS
burn_to_apo = ts_escape.Segments.Insert('eVASegmentTypeManeuver','Burn to Apoapsis','-');
% Change to Finite and leave attitude control as velocity vector
burn_to_apo.SetManeuverType('eVAManeuverTypeFinite');
% Change to correct engine and configure thrust efficiency
burn_to_apo.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_to_apo.Maneuver.ThrustEfficiency = 0.213;
burn_to_apo.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and stopping conditions - Use semimajor axis as a stopping condition
burn_propagator = burn_to_apo.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Heliocentric';
duration_stop = burn_propagator.StoppingConditions.Item('Duration');
duration_stop.Properties.Trip = 1.577e7;
burn_to_apo.Results.Add('Maneuver/DeltaV');
% Set control parameters for target sequence
burn_to_apo.EnableControlParameter('eVAControlManeuverFiniteThrustEfficiency');

%% PROPAGATE - Coast
coast = ts_escape.Segments.Insert('eVASegmentTypePropagate','Coast','-');
% Configure propagator and stopping condition
coast.PropagatorName = 'Heliocentric';
duration_stop = coast.StoppingConditions.Item('Duration');
duration_stop.Properties.Trip = 6.311e6;

%% FINITE MANEUVER TO EARTH ESCAPE
burn_to_escape = ts_escape.Segments.Insert('eVASegmentTypeManeuver','Burn to Earth Escape','-');
% Change to Finite and set attitude thrust vector
burn_to_escape.SetManeuverType('eVAManeuverTypeFinite');
burn_to_escape.Maneuver.SetAttitudeControlType('eVAAttitudeControlTimeVarying');
burn_to_escape.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite VNC(Sun)';
burn_to_escape.Maneuver.AttitudeControl.Az0 = 58.85;
burn_to_escape.Maneuver.AttitudeControl.El0 = -5.162;
% Change to correct engine and configure thrust efficiency
burn_to_escape.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_to_escape.Maneuver.ThrustEfficiency = 0.237;
burn_to_escape.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and stopping conditions - Use semimajor axis as a stopping condition
burn_propagator = burn_to_escape.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Heliocentric';
burn_propagator.MaxPropagationTime = 86400*300;
duration_stop = burn_propagator.StoppingConditions.Item('Duration');
duration_stop.Properties.Trip = 5.004e6;
burn_to_escape.Results.Add('Maneuver/DeltaV');
% Set results to use in target sequence
burn_to_escape.Results.Add('MultiBody/Delta Declination');
burn_to_escape.Results.Item('Delta_Declination').CentralBodyName = 'Earth';
burn_to_escape.Results.Add('MultiBody/Delta Right Asc');
burn_to_escape.Results.Item('Delta_Right_Asc').CentralBodyName = 'Earth';
burn_to_escape.Results.Add('Spherical Elems/R Mag');
burn_to_escape.Results.Add('Target Vector/C3 Energy');
% Set control parameters for target sequence
burn_to_escape.EnableControlParameter('eVAControlManeuverFiniteAz0');
burn_to_escape.EnableControlParameter('eVAControlManeuverFiniteEl0');
burn_to_escape.EnableControlParameter('eVAControlManeuverFiniteThrustEfficiency');
duration_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');

%% Configure TS - Target Earth Escape
% Get handle to differential corrector used in target sequence
dc = ts_escape.Profiles.Item('Differential Corrector');
% Configure Thrust Efficiency control parameters - change their max steps & perturbations
burn1eff_control = dc.ControlParameters.GetControlByPaths('Burn_to_Apoapsis','FiniteMnvr.Thrusting.ThrustEfficiency');
burn1eff_control.Enable = true;
burn1eff_control.MaxStep = 0.01;
burn1eff_control.Perturbation = 0.0001;
burn2eff_control = dc.ControlParameters.GetControlByPaths('Burn_to_Earth_Escape','FiniteMnvr.Thrusting.ThrustEfficiency');
burn2eff_control.Enable = true;
burn2eff_control.MaxStep = 0.005;
burn2eff_control.Perturbation = 0.0001;
% Configure Duration control parameter - change max step and perturbation
sma_control = dc.ControlParameters.GetControlByPaths('Burn_to_Earth_Escape','FiniteMnvr.StoppingConditions.Duration.TripValue');
sma_control.Enable = true;
sma_control.MaxStep = 36000;
sma_control.Perturbation = 60;
% Enable impulsive magnitude control (no need to configure)
dc.ControlParameters.GetControlByPaths('Impulsive_at_Periapsis','ImpulsiveMnvr.Pointing.Spherical.Magnitude').Enable = true;
% Enable thrust vector controls (no need to configure)
pointing_string = 'FiniteMnvr.Pointing.TimeVarying';
dc.ControlParameters.GetControlByPaths('Burn_to_Earth_Escape',sprintf('%s.Az0',pointing_string)).Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_to_Earth_Escape',sprintf('%s.El0',pointing_string)).Enable = true;
% Enable results and configure tolerances
delta_dec_result = dc.Results.GetResultByPaths('Burn_to_Earth_Escape','Delta_Declination');
delta_dec_result.Enable = true;
delta_ra_result = dc.Results.GetResultByPaths('Burn_to_Earth_Escape','Delta_Right_Asc');
delta_ra_result.Enable = true;
rmag_result = dc.Results.GetResultByPaths('Burn_to_Earth_Escape','R_Mag');
rmag_result.Enable = true;
rmag_result.Tolerance = 1e6;
c3_result = dc.Results.GetResultByPaths('Burn_to_Earth_Escape','C3_Energy');
c3_result.Enable = true;
c3_result.Tolerance = 2;
% Set final DC and targeter properties and run modes
dc.MaxIterations = 50;
dc.EnableDisplayStatus = true;
dc.Mode = 'eVAProfileModeIterate';

%% Rearrange MCS to put Backwards Sequence first
MCS.InsertCopy(MCS.Item(5),"Target Asteroid Plane");
MCS.Remove("Backwards from EGA");
MCS.Item(0).Name = "Backwards from EGA";

%% Call function to systematically achieve convergence with all target sequences

fprintf('Trajectory from EGA to asteroid converged! Running Earth escape to EGA.\n\n')

[sat,diverged] = run_target_sequences(sat,"EGA");
















end