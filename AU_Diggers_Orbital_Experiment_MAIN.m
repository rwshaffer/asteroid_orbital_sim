%% Documentation and clear/clc
%AU Digger Orbital Simulation Experiment MAIN
%02/05/2022 Alexander Petsopoulos
%Much of this script is adapted from AGI's
%'STK_Matlab_Object_Test.m'
clear, clc

%% 3. Initialize STK instance/scenario
%grab an already open instance of STK
uiapp = actxGetRunningServer('STK12.Application');

%turn off graphics

% Attach to the STK Object Model
root = uiapp.Personality2;

%%  4. Populate simulation objects: Earth, asteroid, electric propulsion engine

%going to try and take care of this in STK beforehand
%https://agiweb.secure.force.com/faqs/articles/HowTo/interplanetary-ephemeris-summary

%% 5. Insert satellite and configure some settings common to all simulations
%- Turn off any unnecessary graphics to optimize sim. speed

satName = 'BaselineTraj_SemiFinite'; %can change iteratively if necessary
sat = root.CurrentScenario.Children.Item(satName); %satellite handle
ASTK = sat.Propagator; %propagator handle
MCS = ASTK.MainSequence; %Create a handle to the MCS