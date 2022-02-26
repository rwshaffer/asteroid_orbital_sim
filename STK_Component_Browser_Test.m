% Test using the Object Model to change items in the component browser, such as 
% spacecraft engine models.
%
% Some code stolen from STK_Matlab_Object_Test.m
%
% Information on accessing the Component Browser (where engine models and central bodies
% are defined in STK) here:
% https://help.agi.com/stkdevkit/index.htm#stkGator/AstrogatorIntroObjectModel.htm?TocPath=Using%2520Core%2520Libraries%257CSTK%2520Object%2520Model%257CSTK%2520Astrogator%257C_____2
% https://help.agi.com/stkdevkit/Content/DocX/STKObjects~IAgScenario.html
% https://help.agi.com/stkdevkit/Content/DocX/STKObjects~IAgComponentInfoCollection.html
% 


%% Set up STK, insert Astrogator satellite (stolen from STK_Matlab_Object_Test.m)
% Launch a new instance of STK
uiapp = actxserver('STK12.Application');
% or to grab an already open instance of STK
%uiapp = actxGetRunningServer('STK11.Application');

% Attach to the STK Object Model
root = uiapp.Personality2;

% Create a new scenario
root.NewScenario('Components_Test')
% if one already exists, you can access it as in the examples below

% Create a new satellite. See STK Programming Interface Help to see that
% the enumeration for a Satellite object is 'eSatellite' with a value of 18
sat = root.CurrentScenario.Children.New(18, 'ASTG_Sat');
% or connect to an already existing satellite
%sat = root.CurrentScenario.Children.Item('Satellite1');

% Set the new Satellite to use Astrogator as the propagator
sat.SetPropagatorType('ePropagatorAstrogator')
% Note that Astrogator satellites by default start with one Initial State
% and one Propagate segment

% Create a handle to the Astrogator portion of the satellites object model
% for convenience
ASTG = sat.Propagator;


%% Connect to component browser and folder with engine models
comp_brows = root.CurrentScenario.ComponentDirectory.GetComponents('eComponentAstrogator');
engine_models = comp_brows.GetFolder('Engine Models');

%% Clone "Constant Thrust and Isp" engine and change properties
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
ML.sma = 1.9033e+08;
ML.ecc = 0.136364;
ML.inc = 4.37803;
ML.raan = 104.33806;
ML.aop = 183.360698;
ML.M = 240.65329;
ML.period = 524.15224 * 86400; % Units: seconds

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


%% Insert Earth and 1989 ML as planets
earth_obj = root.CurrentScenario.Children.New('ePlanet', 'Earth');
earth_obj.PositionSourceData.set('CentralBody','Earth');
ML_obj = root.CurrentScenario.Children.New('ePlanet', '1989_ML');
ML_obj.PositionSourceData.set('CentralBody','1989_ML');

