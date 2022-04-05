function [sat,diverged] = AU_Diggers_Orbital_Experiment_R2E(sat,mining_duration)

%% CLEARING PARAMETERS
ASTG = sat.Propagator;

ASTG.Options.DrawTrajectoryIn3D = 0; % Turn off drawing trajectory while calculating

MCS = ASTG.MainSequence;
% If this isn't the first return trajectory, we must remove the return target sequence already in the MCS
try
    MCS.Remove('Target Return to Earth')
end

%% TARGET SEQUENCE: Return to Earth
ts_R2E = MCS.Insert('eVASegmentTypeTargetSequence','Target Return to Earth','-');
%% PROPAGATE: Proximity Operations
prox_ops = ts_R2E.Segments.Insert('eVASegmentTypePropagate','Proximity Operations','-');
prox_ops.PropagatorName = 'Heliocentric';
duration_stop = prox_ops.StoppingConditions.Item('Duration');
duration_stop.Properties.Trip = mining_duration * 86400; % Input, mining_duration, converted to seconds
% Set control parameters for target sequence
duration_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');

%% MANEUVER: Burn to Leave Asteroid
burn_leave_ast = ts_R2E.Segments.Insert('eVASegmentTypeManeuver','Burn_to_Leave_Asteroid','-');
% Change to Finite and set attitude to opposite the velocity vector
burn_leave_ast.SetManeuverType('eVAManeuverTypeFinite');
burn_leave_ast.Maneuver.SetAttitudeControlType('eVAAttitudeControlAntiVelocityVector');
% Change to correct engine and configure thrust efficiency
burn_leave_ast.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_leave_ast.Maneuver.ThrustEfficiency = 1;
burn_leave_ast.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and duration stopping condition
burn_propagator = burn_leave_ast.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Heliocentric';
duration_stop = burn_propagator.StoppingConditions.Item('Duration');
duration_stop.Properties.Trip = 2.672e6;
% Set control parameters for target sequence
burn_leave_ast.EnableControlParameter('eVAControlManeuverFiniteThrustEfficiency');
duration_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');

%% PROPAGATE: Coast
prop_coast = ts_R2E.Segments.Insert('eVASegmentTypePropagate','Coast','-');
prop_coast.PropagatorName = 'Heliocentric';
duration_stop = prop_coast.StoppingConditions.Item('Duration');
duration_stop.Properties.Trip = 3.614e6;
% Set control parameters for target sequence
duration_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');

%% MANEUVER: Burn to Earth Return
burn_earth_return = ts_R2E.Segments.Insert('eVASegmentTypeManeuver','Burn_to_Earth_Return','-');
% Change to Finite and set attitude thrust vector
burn_earth_return.SetManeuverType('eVAManeuverTypeFinite');
burn_earth_return.Maneuver.SetAttitudeControlType('eVAAttitudeControlTimeVarying');
burn_earth_return.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite VNC(Sun)';
burn_earth_return.Maneuver.AttitudeControl.Az0 = 255;
burn_earth_return.Maneuver.AttitudeControl.El0 = -6.793;
% Change to correct engine and configure thrust efficiency
burn_earth_return.Maneuver.SetPropulsionMethod('eVAPropulsionMethodEngineModel','AU Diggers Engine');
burn_earth_return.Maneuver.ThrustEfficiency = 0.204;
burn_earth_return.Maneuver.ThrustEfficiencyMode = 'eVAThrustTypeAffectsAccelandMassFlow';
% Configure propagator and Epoch stopping condition
burn_propagator = burn_earth_return.Maneuver.Propagator;
burn_propagator.PropagatorName = 'Heliocentric';
burn_propagator.StoppingConditions.Add('Epoch');
burn_propagator.StoppingConditions.Remove('Duration');
epoch_stop = burn_propagator.StoppingConditions.Item('Epoch');
epoch_stop.Properties.Trip = "11 Jul 2029 13:00:00.00";
% Set control parameters for target sequence
burn_earth_return.EnableControlParameter('eVAControlManeuverFiniteAz0');
burn_earth_return.EnableControlParameter('eVAControlManeuverFiniteEl0');
burn_earth_return.EnableControlParameter('eVAControlManeuverFiniteThrustEfficiency');
epoch_stop.EnableControlParameter('eVAControlStoppingConditionTripValue');
% Configure results for target sequence - Return to Earth
burn_earth_return.Results.Add('MultiBody/Delta Declination');
burn_earth_return.Results.Item('Delta_Declination').CentralBodyName = 'Earth';
burn_earth_return.Results.Add('MultiBody/Delta Right Asc');
burn_earth_return.Results.Item('Delta_Right_Asc').CentralBodyName = 'Earth';
burn_earth_return.Results.Add('Spherical Elems/R Mag');
burn_earth_return.Results.Add('Spherical Elems/V Mag');

%% SETUP FOR TARGET SEQUENCE

% Get handle to differential corrector used in target sequence
dc = ts_R2E.Profiles.Item('Differential Corrector');

% Configure duration control parameters - change max steps
proxops_duration_control = dc.ControlParameters.GetControlByPaths('Proximity_Operations','StoppingConditions.Duration.TripValue');
%proxops_duration_control.Enable = true; % This is an experiment independent variable, so shouldn't be varied in the TS
proxops_duration_control.MaxStep = 86400;
proxops_duration_control.Perturbation = 60;

burn1_duration_control = dc.ControlParameters.GetControlByPaths('Burn_to_Leave_Asteroid','FiniteMnvr.StoppingConditions.Duration.TripValue');
burn1_duration_control.Enable = true;
burn1_duration_control.MaxStep = 86400;
burn1_duration_control.Perturbation = 60;

coast_duration_control = dc.ControlParameters.GetControlByPaths('Coast','StoppingConditions.Duration.TripValue');
coast_duration_control.Enable = true;
coast_duration_control.MaxStep = 86400;
coast_duration_control.Perturbation = 60;

final_epoch_control = dc.ControlParameters.GetControlByPaths('Burn_to_Earth_Return','FiniteMnvr.StoppingConditions.Epoch.TripValue');
final_epoch_control.Enable = true;
final_epoch_control.MaxStep = 86400;
final_epoch_control.Perturbation = 60;

% Configure thrust efficiency control parameter - change max step & perturbation
burn2_eff_control = dc.ControlParameters.GetControlByPaths('Burn_to_Earth_Return','FiniteMnvr.Thrusting.ThrustEfficiency');
burn2_eff_control.Enable = true;
burn2_eff_control.MaxStep = 0.005;
burn2_eff_control.Perturbation = 0.0001;

% Configure maneuver pointing controls - no need to change steps
dc.ControlParameters.GetControlByPaths('Burn_to_Earth_Return','FiniteMnvr.Pointing.TimeVarying.Az0').Enable = true;
dc.ControlParameters.GetControlByPaths('Burn_to_Earth_Return','FiniteMnvr.Pointing.TimeVarying.El0').Enable = true;


% Enable results and configure tolerances
delta_dec_result = dc.Results.GetResultByPaths('Burn_to_Earth_Return','Delta_Declination');
delta_dec_result.Enable = true;
delta_dec_result.Tolerance = 1;
delta_ra_result = dc.Results.GetResultByPaths('Burn_to_Earth_Return','Delta_Right_Asc');
delta_ra_result.Enable = true;
delta_ra_result.Tolerance = 0.5;
rmag_result = dc.Results.GetResultByPaths('Burn_to_Earth_Return','R_Mag');
rmag_result.Enable = true;
rmag_result.Tolerance = 2e6;
vmag_result = dc.Results.GetResultByPaths('Burn_to_Earth_Return','V_Mag');
%vmag_result.Enable = true; % Better when this isn't used
vmag_result.Tolerance = 3;

% Set final DC and targeter properties and run modes
dc.MaxIterations = 50;
dc.EnableDisplayStatus = true;
dc.Mode = 'eVAProfileModeIterate';

[sat,diverged] = run_target_sequences(sat);

end
