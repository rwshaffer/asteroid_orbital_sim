%% Documentation and clear/clc
%AU Digger Orbital Simulation Experiment MAIN
%02/05/2022 Alexander Petsopoulos
%Much of this script is adapted from AGI's
%'STK_Matlab_Object_Test.m'
clear, clc

%% 3. Initialize STK instance/scenario
% Launch a new instance of STK
uiapp = actxserver('STK12.Application');

% or to grab an already open instance of STK
%uiapp = actxGetRunningServer('STK11.Application');

%turn off graphics

% Attach to the STK Object Model
root = uiapp.Personality2;

% Create a new scenario
root.NewScenario('OrbitalSim')
% if one already exists, you can access it as in the examples below

%%  4. Populate simulation objects: Earth, asteroid, electric propulsion engine

%Initializing Jupiter as a test for how to create a planet
% IAgScenario scenario: Scenario object
planet = root.CurrentScenario.Children.New('ePlanet', 'Jupiter');
planet.CommonTasks.SetPositionSourceCentralBody('Jupiter', 'eEphemJPLDE');

% IAgPlanet planet: Planet object
planet2D = planet.Graphics;
planet2D.Color = 255;   % Red
planet2D.Inherit = false;
planet2D.OrbitVisible = true;
planet2D.SubPlanetPointVisible = false;
planet2D.SubPlanetLabelVisible = false;

%% 5. Insert satellite and configure some settings common to all simulations

%- Turn off any unnecessary graphics to optimize sim. speed

% Create our satellite.
sat = root.CurrentScenario.Children.New('eSatellite', 'DiggerSat');
% or connect to an already existing satellite
%sat = root.CurrentScenario.Children.Item('DiggerSat');

% Set the new Satellite to use Astrogator as the propagator
sat.SetPropagatorType('ePropagatorAstrogator')
% Note that Astrogator satellites by default start with one Initial State
% and one Propagate segment

% Create a handle to the Astrogator portion of the satellites object model
% for convenience
DigPropagator = sat.Propagator;

% In MATLAB, you can use the .get command to return a list of all
% "attributes" or properties of a given object class. Examine the
% Astrogator Object Model Diagram to see a depiction of these.
DigPropagator.get
%    MainSequence: [1x1 Interface.AGI_STK_Astrogator_9.IAgVAMCSSegmentCollection]
%         Options: [1x1 Interface.AGI_STK_Astrogator_9._IAgVAMCSOptions]
%    AutoSequence: [1x1 Interface.AGI_STK_Astrogator_9.IAgVAAutomaticSequenceCollection]

% In MATLAB, you can use the .invoke command to return a list of all
% "methods" or functions of a given object class. Examine the Astrogator
% Object Model Diagram to see a depiction of these.
DigPropagator.invoke
%     RunMCS = void RunMCS(handle)
%     BeginRun = void BeginRun(handle)
%     EndRun = void EndRun(handle)
%     ClearDWCGraphics = void ClearDWCGraphics(handle)
%     ResetAllProfiles = void ResetAllProfiles(handle)
%     ApplyAllProfileChanges = void ApplyAllProfileChanges(handle)
%     AppendRun = void AppendRun(handle)
%     AppendRunFromTime = void AppendRunFromTime(handle, Variant, AgEVAClearEphemerisDirection)
%     AppendRunFromState = void AppendRunFromState(handle, handle, AgEVAClearEphemerisDirection)
%     RunMCS2 = AgEVARunCode RunMCS2(handle)

% At any place in the STK or Astrogator OM, use the .get or .invoke
% commands to inspect the structure of the object model and help find the
% desired properties or methods

%%%
% Adding and Removing segments
%%%

% Collections
% In the OM, groupings of the same kind of object are referred to as
% Collections. Examples include Sequences (including the MainSequence and
% Target Sequences) which hold groups of segments, Segments which may hold
% groups of Results, and Propagate Segments which may hold groups of
% Stopping Conditions.
% In general, all Collections have some similar properties and methods and
% will be interacted with the same way. The most common elements of a
% Collection interface are
%   Item(argument) - returns a handle to a particular element of
%   the collection
%   Count - the number of elements in this collection
%   Add(argument) or Insert(argument) - adds new elements to the collection
%   Remove, RemoveAll - removes elements from the collection
% Other methods like Cut, Copy, and Paste may be available depending on the
% kind of collection

% Create a handle to the MCS and remove all existing segments
MCS = DigPropagator.MainSequence;
MCS.RemoveAll;

% Functions can also be called directly without needing to create a
% separate handle. This will also work:
% ASTG.MainSequence.RemoveAll;

%%% Define the Initial State %%%

% Use the Insert method to add a new Initial State to the MCS. The Insert
% method requires an enumeration as one of its arguments. Enumerations are
% a set of pre-defined options for certain methods and can be found in the
% Help for that given method.
MCS.Insert('eVASegmentTypeInitialState','Inner Orbit','-');

% The Insert command will also return a handle to the segment it creates
propagate = MCS.Insert('eVASegmentTypePropagate','Propagate','-');

%%%
% Configuring Segment properties
%%%

% Create a handle to the Initial State Segment, set it to use Modified
% Keplerian elements and assign new initial values
initstate = MCS.Item('Inner Orbit');
initstate.SetElementType('eVAElementTypeTargetVectorOutgoingAsymptote');