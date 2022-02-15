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

satName = 'BaselineTraj_SemiFinite';
sat = root.CurrentScenario.Children.Item(satName);
ASTK = sat.Propagator;

% Create a handle to the MCS and remove all existing segments
MCS = ASTK.MainSequence;
MCS.RemoveAll;
MCS.Insert('eVASegmentTypeInitialState','Inner Orbit','-');

% The Insert command will also return a handle to the segment it creates
propagate = MCS.Insert('eVASegmentTypePropagate','Propagate','-');

% Create a handle to the Initial State Segment, set it to use Modified
% Keplerian elements and assign new initial values
initstate = MCS.Item('Inner Orbit');
initstate.SetElementType('eVAElementTypeTargetVectorOutgoingAsymptote');