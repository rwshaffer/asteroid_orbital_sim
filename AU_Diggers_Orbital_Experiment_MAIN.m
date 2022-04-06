%% Documentation and clear/clc
%AU Digger Orbital Simulation Experiment MAIN
%Much of this script is adapted from AGI's
%'STK_Matlab_Object_Test.m'

% Information on accessing the Component Browser (where engine models and central bodies
% are defined in STK) here:
% https://help.agi.com/stkdevkit/index.htm#stkGator/AstrogatorIntroObjectModel.htm?TocPath=Using%2520Core%2520Libraries%257CSTK%2520Object%2520Model%257CSTK%2520Astrogator%257C_____2
% https://help.agi.com/stkdevkit/Content/DocX/STKObjects~IAgScenario.html
% https://help.agi.com/stkdevkit/Content/DocX/STKObjects~IAgComponentInfoCollection.html
%

% inputVec will be the vector of inputs, taking data from Terelle's GUI
clear, clc



%% User Inputs

launch_dates = ["7 Jul 2025 07:00:00.000",...
                "16 Jul 2025 07:00:00.000"];

mining_durations = [482,500]; % Days


outbound_traj = 'Direct'; % Options are 'EGA' or 'Direct'

asteroid_uncertainty = 0; % Percent orbit uncertainty to apply to all elements uniformly




%% 3. Initialize STK instance/scenario

uiapp = actxserver('STK12.Application');

%turn off graphics

% Attach to the STK Object Model
root = uiapp.Personality2;

%Create scenario
root.NewScenario('Orbital_Sim')
%%  4. Populate simulation objects: Earth, asteroid, electric propulsion engine

%% Connect to component browser and folder with engine models
comp_brows = root.CurrentScenario.ComponentDirectory.GetComponents('eComponentAstrogator');
engine_models = comp_brows.GetFolder('Engine Models');

%Clone "Constant Thrust and Isp" engine and change properties
engine_models.DuplicateComponent('Constant Thrust and Isp','AU Diggers Engine');
our_engine = engine_models.Item('AU Diggers Engine');

thrust = .101; %Newtons
Isp = 1710; % seconds

our_engine.set('Thrust',thrust);
our_engine.set('Isp',Isp);

%% Connect to folder with planets
central_bodies = comp_brows.GetFolder('Central Bodies');

%% Create 1989 ML based on cloning asteroid template
central_bodies.DuplicateComponent('AsteroidTemplate','1989 ML');
ML_body = central_bodies.Item('1989 ML');

%% Find 1989 ML orbital parameters from JPL
% Note: Probably easiest to save these as a .mat file
ML.epoch = 2.4596e+06;
ML.sma = 1.9033e+08 * (1+asteroid_uncertainty/100);
ML.ecc = 0.136364 * (1+asteroid_uncertainty/100);
ML.inc = 4.37803 * (1+asteroid_uncertainty/100);
ML.raan = 104.33806 * (1+asteroid_uncertainty/100);
ML.aop = 183.360698 * (1+asteroid_uncertainty/100);
ML.M = 240.65329 * (1+asteroid_uncertainty/100);
ML.period = 524.15224 * 86400 * (1+asteroid_uncertainty/100); % Units: seconds

% Define a few other quantities used by STK to define the orbit
ML.longofperiapsis = ML.raan + ML.aop; % Longitude of periapsis (deg)
ML.meanlongitude = ML.raan + ML.aop + ML.M; % Mean longitude (deg)
ML.meanlongrate = 360/ML.period; % Mean angular rate (deg/sec)

%% Define new asteroid object's orbital parameters
ML_orbit = ML_body.DefaultEphemerisData;
%ML_orbit.get % View orbital parameter variable names and values
ML_orbit.Epoch = ML.epoch;
ML_orbit.SemiMajorAxis = ML.sma;
ML_orbit.Eccentricity = ML.ecc;
ML_orbit.Inclination = ML.inc;
ML_orbit.RAAN = ML.raan;
ML_orbit.ArgOfPeriapsis = ML.longofperiapsis; % NOTE: This variable is called ArgOfPeriapsis, but it seems to refer
%                                               to longitude of periapsis, which is different. Be careful!
ML_orbit.MeanLongitude = ML.meanlongitude;
ML_orbit.MeanLongitudeRate = ML.meanlongrate;

ML_body.SetDefaultEphemerisByName('Analytic Orbit');

%% Insert Earth and 1989 ML as planets
earth_obj = root.CurrentScenario.Children.New('ePlanet', 'Earth');
earth_obj.PositionSourceData.set('CentralBody','Earth');
ML_obj = root.CurrentScenario.Children.New('ePlanet', '1989_ML');
ML_obj.PositionSourceData.set('CentralBody','1989_ML');



%% 5. Insert satellite and configure some settings common to all simulations
%- Turn off any unnecessary graphics to optimize sim. speed


% Create a new satellite. See STK Programming Interface Help to see that
% the enumeration for a Satellite object is 'eSatellite' with a value of 18
sat = root.CurrentScenario.Children.New(18, 'AU_Digger_Sat');

% Set the new Satellite to use Astrogator as the propagator
sat.SetPropagatorType('ePropagatorAstrogator')
% Note that Astrogator satellites by default start with one Initial State
% and one Propagate segment

% Create a handle to the Astrogator portion of the satellites object model
% for convenience
ASTG = sat.Propagator;
MCS = ASTG.MainSequence; %Create a handle to the MCS

% Create new axes to track thrust vector with respect to 1989 ML
% VNC (velocity, normal, conormal) coordinate system centered at 1989 ML
axesFactory = sat.Vgt.Axes.Factory;
vncML = axesFactory.Create('VNC(1989_ML)','','eCrdnAxesTypeTrajectory');
vncML.ReferenceSystem.SetPath('CentralBody/1989_ML J2000');
vncML.TrajectoryAxesType = 'eCrdnTrajectoryAxesVVLH';


%% --------------------------------------------------------------------- %%
% Set up results table
numTrajectories = length(launch_dates) * length(mining_durations);
numOutputs = 7; % Hard-coded, based on extract_sat.m
dataTable = zeros(numOutputs,numTrajectories);

traj_num = 0; % Keep track of all trajectories (total number of full trajectories modeled)

% Keep track of inputs given to each trajectory
%traj_type_table = [];
uncertainty_table = zeros(1,numTrajectories);
%launch_date_table = zeros(1,numTrajectories); 
mining_dur_table = zeros(1,numTrajectories);

% LOOP THROUGH ALL USER INPUTS
for i = 1:length(launch_dates)
    %% 6a. Build Mission Control Sequence for satellite outbound traj. based on user inputs
    launch_date = launch_dates(i);
    
    % Display info on outbound trajectory about to be modeled
    fprintf("\n\n--------------------------------------------------------\n")
    fprintf("Outbound Trajectory %d: Launch Date = %s\n\n",i,launch_date)
    
    switch outbound_traj
        case 'Direct'
            [sat,diverged] = escape_to_earth_plane_MCS(ML,sat,launch_date);
        case 'EGA'
            fprintf("EGA modeling is incomplete, try again later.\n")
    end
    %% 6b. Assuming outbound traj. converged, build in the return trajectory
    if diverged
        fprintf('Outbound trajectory diverged! Moving to next outbound traj.:\n\n')
    else
        fprintf('Outbound trajectory converged! Building return trajectory(s):\n\n')
    end
    
    for j = 1:length(mining_durations)
        % Update new mining duration and # trajectory modeled
        traj_num = traj_num + 1;
        mining_duration = mining_durations(j);

        % Save inputs to this particular trajectory (for final table)
        traj_type_table(traj_num) = string(outbound_traj);
        uncertainty_table(traj_num) = asteroid_uncertainty;
        launch_date_table(traj_num) = string(launch_date);
        mining_dur_table(traj_num) = mining_duration;
        
        if ~diverged
            % Display info on return trjaectory about to be modeled
            fprintf("-----------------------------------------\n")
            fprintf("Return Trajectory %d: Mining Duration = %d days\n\n",j,mining_duration)
            [sat,diverged] = AU_Diggers_Orbital_Experiment_R2E(sat,mining_duration);
            if ~diverged
                %% 7. Extracting data
                fprintf("Return trajectory converged! SUMMARY: \n")
                dataTable(:,traj_num) = extract_data(sat,mining_duration)';
            else
                dataTable(:,traj_num) = nan;
            end
        else
            dataTable(:,traj_num) = nan;
        end
    
    end
  
end

%% Saving Data
fprintf('All trajectories are complete! Final results:\n\n')

% Inputs for each trajectory
InputLaunchDates = launch_date_table';
InputMiningDurations = mining_dur_table';
InputTrajType = traj_type_table';
InputUncertainty = uncertainty_table';
% Outputs for each trajectory
OutboundFinDV = dataTable(1,:)';
OutboundImpDV = dataTable(2,:)';
ReturnDV = dataTable(3,:)';
MaxThrust = dataTable(4,:)';
C3Energy = dataTable(5,:)';
OutboundTimeOfFlight = dataTable(6,:)';
ReturnTimeOfFlight = dataTable(7,:)';
% Consolidate inputs and outputs into a table and display it
outputTable = table(InputLaunchDates,InputMiningDurations,InputTrajType,InputUncertainty,...
    OutboundFinDV,OutboundImpDV,ReturnDV,MaxThrust,C3Energy,OutboundTimeOfFlight,ReturnTimeOfFlight);
disp(outputTable)
% Write outputs to Excel
%writetable(outputTable);
filename = 'Output Data.xlsx';
save_data = input('Enter "y" to save data. Warning: this deletes all other entries in file!\n',"s");
if save_data == 'y'
    writetable(outputTable,filename,'WriteMode','replacefile');
    fprintf('Saved!\n')
end









