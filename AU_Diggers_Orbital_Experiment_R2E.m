%% CLEARING PARAMETERS
ASTG = sat.Propagator;

ASTG.Options.DrawTrajectoryIn3D = 0; % Turn off drawing trajectory while calculating

MCS = ASTG.MainSequence;
MCS.RemoveAll;
%% TARGET SEQUENCE: Return to Earth
ts_R2E = MCS.Insert('eVASegmentTypeTargetSequence','Return to Earth','-');
%% INITAL STATE: Setting Initial State

%Parameters InitialState/Elements
initstate = ts_R2E.Segments.Insert('eVASegmentTypeInitialState','InitialState','-');
initstate.CoordSystemName = "CentralBody/1989_ML J2000";
initstate.OrbitEpoch = "1 Jan 2022 17:00:00.000";
initstate.Element.X = 0;
initstate.Element.Y = 0;
initstate.Element.Z = 0;
initstate.Element.Vx = 0;
initstate.Element.Vy = 0;
initstate.Element.Vz = 0;

% Set initial S/C parameters
initstate.SpacecraftParameters.DryMass = 224;
initstate.FuelTank.FuelMass = 170;

% Set control parameters for target sequence
initstate.EnableControlParameter('eVAControlInitStateEpoch')
%% MANEUVER: Time Varying Burn

Maneuver1 = ts_R2E.Segments.Insert('eVASegmentTypeManeuver','Maneuver','-');
Maneuver1.SetManeuverType('eVAManeuverTypeFinite');
Maneuver1.Maneuver.Propagator.PropagatorName = 'Heliocentric';
Maneuver1.Maneuver.SetAttitudeControlType('eVAAttitudeControlTimeVarying');
Maneuver1.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite VNC(Sun)';
Maneuver1.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');

Maneuver1.Maneuver.Propagator.StoppingConditions.Add('Epoch');
Maneuver1.Maneuver.Propagator.StoppingConditions.Remove('Duration');
sma_stop = Maneuver1.Maneuver.Propagator.StoppingConditions.Item('Epoch');
sma_stop.Properties.Trip = "1 Jul 2022 17:00:00.00";

Maneuver1.Maneuver.AttitudeControl.Az0 = 180;

% Set control parameters for target sequence
Maneuver1.EnableControlParameter('eVAControlManeuverFiniteAz0');
Maneuver1.EnableControlParameter('eVAControlManeuverFiniteEl0');
sma_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');

% Configure results for target sequence - Return to Earth
Maneuver1.Results.Add('MultiBody/Delta Declination');
Maneuver1.Results.Item('Delta_Declination').CentralBodyName = 'Earth';
Maneuver1.Results.Add('MultiBody/Delta Right Asc');
Maneuver1.Results.Item('Delta_Right_Asc').CentralBodyName = 'Earth';
Maneuver1.Results.Add('Spherical Elems/R Mag');
Maneuver1.Results.Item('R_Mag').ReferencePointName = 'CentralBody/1989_ML Center';

%% SETUP FOR TARGET SEQUENCE

% Get handle to differential corrector used in target sequence
dc = ts_R2E.Profiles.Item('Differential Corrector');

% Configure Epoch control parameter - change its max step & perturbation
init_epoch_control = dc.ControlParameters.GetControlByPaths('InitialState','InitialState.Epoch');
init_epoch_control.Enable = true;
init_epoch_control.MaxStep = 20000;
init_epoch_control.Perturbation = 500;

fin_epoch_control = dc.ControlParameters.GetControlByPaths('Maneuver','FiniteMnvr.StoppingConditions.Epoch.TripValue');
fin_epoch_control.Enable = true;
fin_epoch_control.MaxStep = 20000;
fin_epoch_control.Perturbation = 500;

dc.ControlParameters.GetControlByPaths('Maneuver','FiniteMnvr.Pointing.TimeVarying.Az0').Enable = true;
dc.ControlParameters.GetControlByPaths('Maneuver','FiniteMnvr.Pointing.TimeVarying.El0').Enable = true;


% Enable results and configure tolerances
delta_dec_result = dc.Results.GetResultByPaths('Maneuver','Delta_Declination');
delta_dec_result.Enable = true;
delta_dec_result.Tolerance = 0.1;
delta_ra_result = dc.Results.GetResultByPaths('Maneuver','Delta_Right_Asc');
delta_ra_result.Enable = true;
delta_ra_result.Tolerance = 0.1;
rmag_result = dc.Results.GetResultByPaths('Maneuver','R_Mag');
rmag_result.Enable = true;
rmag_result.Tolerance = 1.5e6;

% Set final DC and targeter properties and run modes
dc.MaxIterations = 50;
dc.EnableDisplayStatus = true;
dc.Mode = 'eVAProfileModeIterate';
ts_R2E.Action = 'eVATargetSeqActionRunActiveProfiles';

ASTG.RunMCS

ASTG.RunMCS
