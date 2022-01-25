function [final_orbit,B_dot] = earth_grav_assist(initial_orbit, earth)
% Construct an Earth gravity assist that changes the S/C heliocentric orbit
% inclination to a desired value.
% See info on B-plane targeting here: https://ai-solutions.com/_freeflyeruniversityguide/the_b_plane.htm
%
% INPUTS
% initial_orbit: struct containing S/C heliocentric orbital elements at beginning of gravity assist.
% earth: struct containing Earth's orbital elements at the time of the gravity assist.
%   NOTE: each of the above structs contains the following fields:
%   epoch - epoch at which orbit is based (JD)
%   sma - semimajor axis (km)
%   ecc - eccentricity
%   inc - inclination (rad)
%   raan - right ascension of the ascending node (rad)
%   aop - argument of periapsis (rad)
%   nu - true anomaly (rad)
% desired_inclination: inclination of final heliocentric orbit (rad)
%
% OUTPUTS
% final_orbit: struct (same fields as above) containing S/C heliocentric orbital elements after gravity assist
% B_dot: Target vector [B_dot_T; B_dot_R] on the S/C B-plane during flyby
%

%% Constants
mu_earth = 3.986e5; % Earth's gravitational constant (km^3/s^2)

%% Compute nitial Cartesian state relative to Earth
[sc_helio_pos_init,sc_helio_vel_init] = Kepler2Carts(initial_orbit,'sun'); % S/C final Cartesian state relative to sun
[earth_helio_pos,earth_helio_vel] = Kepler2Carts(earth,'sun'); % Earth's state relative to sun
% Convert S/C state from rel. to sun to rel. to Earth
sc_earth_pos_init = sc_helio_pos_init - earth_helio_pos; % S/C position relative to Earth
sc_earth_vel_init = sc_helio_vel_init - earth_helio_vel; % S/C velocity relative to Earth
% sc_earth_cart_init.X = sc_helio_cart_init.X - earth_helio_cart.X;
% sc_earth_cart_init.Y = sc_helio_cart_init.Y - earth_helio_cart.Y;
% sc_earth_cart_init.Z = sc_helio_cart_init.Z - earth_helio_cart.Z;
% sc_earth_cart_init.Vx = sc_helio_cart_init.Vx - earth_helio_cart.Vx;
% sc_earth_cart_init.Vy = sc_helio_cart_init.Vy - earth_helio_cart.Vy;
% sc_earth_cart_init.Vz = sc_helio_cart_init.Vz - earth_helio_cart.Vz;

%% Hyperbolic Keplerian elements (rel. to Earth) at beginning of gravity assist
[hyp_elements,e_vec,h_vec,eps] = Carts2Kepler(sc_earth_pos_init,sc_earth_vel_init,'earth');
% Unpack orbital elements into individual variables for clarity
sma = hyp_elements.sma;
ecc = hyp_elements.ecc;
inc = hyp_elements.inc;
raan = hyp_elements.raan;
aop = hyp_elements.aop;
nu = hyp_elements.nu;

%% Determine inbound asymptotic velocity rel. to Earth
Vinf = sqrt(2*eps); % Asymptotic velocity magnitude (km/s)
Vinf_hat_in = sc_earth_vel_init/norm(sc_earth_vel_init); % Direction of initial velocity
%   Note: this approximates the asymptotic velocity as identical to the initial velocity
Vinf_vec_in = Vinf * Vinf_hat_in; % Asymptotic velocity vector (km/s)
%Vinf_vec_in = sc_earth_vel_init;
%Vinf = norm(Vinf_vec_in); % Magnitude of asymptotic velocity (same for outbound and inbound)
%Vinf_hat_in = Vinf_vec_in/Vinf; % Unit vector in direction of incoming asmyptotic velocity

%% Switch between inertial frame (IJK) and orbital frame (PQW)
R_PQW_IJK = [cos(raan)*cos(aop)-sin(raan)*sin(aop)*cos(inc) -cos(raan)*sin(aop)-sin(raan)*cos(aop)*cos(inc) sin(raan)*sin(inc);
    sin(raan)*cos(aop)+cos(raan)*sin(aop)*cos(inc) -sin(raan)*sin(aop)+cos(raan)*cos(aop)*cos(inc) -cos(raan)*sin(inc);
    sin(aop)*sin(inc) cos(aop)*sin(inc) cos(inc)]; % Rotation from PQW to IJK
R_IJK_PQW = R_PQW_IJK'; % Rotation from IJK to PQW


Vinf_in_PQW = R_IJK_PQW * Vinf_vec_in; % Initial velocity in PQW frame
Rinf_in_PQW = R_IJK_PQW * sc_earth_pos_init; % Initial position in PQW frame
% Note: this position vector is labeled "Rinf", but its true position is not negligible

%% Find outbound asymptote in PQw, then IJK frames
% Hyperbolic trajectory is symmetric about e_vec!
Vinf_out_PQW = -Vinf_in_PQW; % Initialize outbound asymptotic velocity
Vinf_out_PQW(2) = -Vinf_out_PQW(2); % Flip sign of conormal (Q) velocity component
Vinf_vec_out = R_PQW_IJK * Vinf_out_PQW; % Oubound velocity in IJK frame

Rinf_out_PQW = Rinf_in_PQW; % Initialize outbound asymptotic position
Rinf_out_PQW(2) = -Rinf_in_PQW(2); % Flip sign of conormal (Q) position component
Rinf_vec_out = R_PQW_IJK * Rinf_out_PQW; % Oubound position in IJK frame


%% Final heliocentric orbit Cartesian state
sc_earth_vel_final = Vinf_vec_out; % S/C final velocity rel. to Earth
sc_earth_pos_final = Rinf_vec_out; % S/C final position rel. to Earth

sc_helio_vel_final = sc_earth_vel_final + earth_helio_vel; % S/C final velocity rel. to sun
sc_helio_pos_final = sc_earth_pos_final + earth_helio_pos; % S/C final position rel. to sun

%% Final heliocentric orbit Keplerian elements
final_orbit = Carts2Kepler(sc_helio_pos_final,sc_helio_vel_final,'sun');

%% Hyperbolic trajectory in Earth's SOI - B-Plane
% Based on this tutorial: https://ai-solutions.com/_freeflyeruniversityguide/the_b_plane.htm

h_hat = h_vec/norm(h_vec); % Angular momentum unit vector
e_hat = e_vec/ecc; % Eccentricity unit vector
beta = acos(1/ecc); % Beta angle (rad)
S_hat = cos(beta)*e_hat + sin(beta) * cross(h_hat,e_hat); % S unit vector (direction of Vinf at beginning of grav. assist)
N_hat = [0;0;1]; % N unit vector (defined as reference frame Z-axis)
T_hat = cross(S_hat,N_hat)./norm(cross(S_hat,N_hat)); % T unit vector (B-plane)
R_hat = cross(S_hat,T_hat); % R unit vector (B-plane)
B_mag = abs(sma) * sqrt(ecc^2 - 1); % Magnitude of S/C B-vector
B_hat = cross(S_hat,h_hat); % Direction of S/C B-vector
B_vec = B_mag * B_hat; % S/C B-vector
B_dot_T = dot(B_vec,T_hat); % B dot T, projection of B vector onto T axis
B_dot_R = dot(B_vec,R_hat); % B dot R, projection of B vector onto R axis
B_dot = [B_dot_T; B_dot_R]; % B dot vector, defining the target point on the B-plane
theta = acos(dot(B_vec,T_hat)./(norm(B_vec) * norm(T_hat))); % Vertex angle, defines the plane of the hyperbolic trajectory
delta = pi - 2*beta; % angle (rad) between incoming and outgoing asymptotic velocities


end