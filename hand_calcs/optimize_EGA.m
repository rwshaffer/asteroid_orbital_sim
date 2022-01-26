function [initial_orbit] = optimize_EGA(initial_orbit,earth,target)
% Approximately optimize the initial orbit required to match a target's orbital parameters
% after an Earth gravity assist (EGA)
%
% Currently, the following orbital parameters are optimized:
%   inclination
%   right ascension of the ascending node
%
% INPUTS
% initial_orbit: struct describing initial, unoptimized orbit
% earth: struct containing orbital elements of Earth at the time of the gravity assist
% target: struct containing orbital elements of target at the time of EGA
%   NOTE: each of the above structs contain the following elements:
%   epoch - epoch at which orbit is based (JD)
%   sma - semimajor axis (km)
%   ecc - eccentricity
%   inc - inclination (rad)
%   raan - right ascension of the ascending node (rad)
%   aop - argument of periapsis (rad)
%   nu - true anomaly (rad)
%
% OUTPUTS
% initial_orbit: struct containing optimized orbital parameters of S/C before EGA that will best
%   match the target's orbital parameters
%

%% Constants/Parameters
deg2rad = pi/180;
rad2deg = 180/pi;

disp_text = false; % Toggle on/off to print out progress of for loops
plots = false; % Toggle on/off to create contour plots


%% Set up input space

% Vary inclination between Earth's inclination (inc. = 0) and target's incliination.
%   Beyond the target's inclination, there is no point conducting a gravity assist,
%   since the delta-V to reach that inclination will be higher than just burning to
%   the target inclination directly.
inclinations = linspace(0,target.inc,200); % Vector containing all inclinations of interest

% Vary RAAN +/- 90 degrees from target's RAAN (there is likely a better way to do this)
raans = linspace(target.raan + pi/2, target.raan + 3*pi/2,360); % Vector containing all RAANs of interest in radians
% NOTE: THIS WORKED FOR 1989 ML BUT NOT EXACTLY SURE WHY. MAY NEED TO LOOK AT ALL 360 DEG. OF RAAN
raans = linspace(0,2*pi,360);


%% Set up nested for loops
% Create matrices used to store final orbital elements (inc. and raan)
final_incs = zeros(length(inclinations),length(raans)); % Matrix of final inclinations
final_raans = zeros(length(inclinations),length(raans)); % Matrix of final RAANs

%% Simulate EGAs for all initial inclinations and RAANs
for i = 1:length(inclinations)
    initial_orbit.inc = inclinations(i); % Update inclination of pre-EGA orbit   
    if disp_text
        fprintf('Round %0.0d, i = %0.2f\n',i,inclinations(i)*rad2deg) % Display current inclination
    end
    for j = 1:length(raans)
        initial_orbit.raan = raans(j); % Update RAAN of pre-EGA orbit
        final_orbit = earth_grav_assist(initial_orbit, earth); % Model gravity assist and find final orbit
        
        final_incs(i,j) = real(final_orbit.inc); % Store final orbit inclination
        final_raans(i,j) = real(final_orbit.raan); % Store final orbit RAAN
    end
end

%% Plot results
if plots
    % Contour plot of final orbit inclination as a function of initial inc. and RAAN
    figure()
    contour_inc_level = [target.inc target.inc] * rad2deg;
    contour(raans*rad2deg,inclinations*rad2deg,final_incs*rad2deg,contour_inc_level,'ShowText','on')
    title('Final Orbit Inclination (deg)')
    xlabel('Initial RAAN (deg)')
    ylabel('Initial Inc. (deg)')

    % Contour plot of final orbit RAAN as a function of initial inc. and RAAN
    figure()
    contour_raan_level = [target.raan target.raan] * rad2deg;
    contour(raans*rad2deg,inclinations*rad2deg,final_raans*rad2deg,contour_raan_level,'ShowText','on')
    title('Final Orbit RAAN (deg)')
    xlabel('Initial RAAN (deg)')
    ylabel('Initial Inc. (deg)')

end

%% Find all orbits that match the target inc./RAAN within a desired tolerance

inc_tol = 1 * deg2rad; % Initial tolerance of diff. between final and target inclination (rad).
raan_tol = 5 * deg2rad; % Initial tolerance of diff. between final and target RAAN (deg)

% Set up while loop to check for valid orbits
orbit_match = zeros(size(final_incs)); % Matrix to track inc/RAAN combos that match target orbit
while nnz(orbit_match) == 0 % Iterate through tolerance checking until a valid orbit is found
    
    inc_difference = abs(final_incs - target.inc); % Matrix of diff. b/w final and target inc. for each initial orbit
    raan_difference = abs(final_raans - target.raan); % Matrix of diff. b/w final and target raan. for each initial orbit

    inc_match = inc_difference < inc_tol; % Elements of inc_difference matrix that meet tolerance
    raan_match = raan_difference < raan_tol; % Elements of raan_difference matrix that meet tolerance
    orbit_match = inc_match .* raan_match; % Elements that meet both tolerances

    % If the number of nonzero elements in orbit_match is 0, 
    % the while loop will continue.
    % Increase the tolerances (for next iteration)
    inc_tol = inc_tol * 1.2;
    raan_tol = raan_tol * 1.2;
end

%% Find the matching orbit with the minimum inclination
% Inclination determines the row of the matrix, with smaller inclinations at the top.
% Iterate through rows of the orbit_match matrix
for i = 1:length(inclinations)
    if any(orbit_match(i,:) == 1) % Check for any 1's in the row, meaning that an orbit at that inclination matches the target
        optimal_inc = inclinations(i); % This initial inclination is optimal
        optimal_raan = raans(find(orbit_match(i,:),1,'first')); % Take lowest initial RAAN that works w/ that inc. to be optimal
        break
    end
end

% Update initial orbit with new optimized parameters
initial_orbit.inc = optimal_inc;
initial_orbit.raan = optimal_raan;

if disp_text
    fprintf('\nInitial inc.: %0.2f deg\n',optimal_inc*rad2deg)
    fprintf('Initial RAAN: %0.2f deg\n',optimal_raan*rad2deg)
end

end