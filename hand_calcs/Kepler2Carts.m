function [position,velocity] = Kepler2Carts(orbit,central_body)
%--------------------------------------------------------------------------------------------------------%
%
% 			USAGE: Conversion of Keplerian Classic Orbital Elements into inertial reference system OXYZ
%
% 			AUTHOR: Thameur Chebbi(PhD)		E-MAIL: chebbythamer@gmail.com
%           Modified by Ryan Shaffer (not PhD)
%
% 			DATE: 01,Oct,2020
%
% 			DESCIPTION:      This function is created to convert the classical 
%               		     orbital elements to cartesian position and velocity 
%		       		         parameters of any satellite orbit in an inertial
%                            reference system centered at the central body and aligned
%                            with the ecliptic plane (sun-centered) or equator (Earth_centered).
%
% 			INPUT:
% 			sma:    Semi-Major Axis..............(Km)							
% 			ecc:    Eccentricity											    
% 			inc:	Inclination..................(rad)							
% 			w:	    Argument of perigee..........(rad)	
% 			nu:	    Satellite position...........(rad)							
% 			RAAN:	Right Asc. of Ascending Node.(rad)							
%
% 			OUTPUT:
%           
%			Position Components: 		
% 			[X; Y; Z]...(Km)
%
%			Velocity Components:
% 			[Vx; Vy; Vz]...(Km/s) 														
%
%
%%---------------------------- Constants ----------------------------------------------%
mu_earth = 3.986 * 10^5; % Earth Gravitational Constant
mu_sun =  1.3271233e11; % Sun Gravitational Constant

%
%%-------------------------- Unpack orbit structure -----------------------------------%
sma = orbit.sma;
ecc = orbit.ecc;
inc = orbit.inc;
raan = orbit.raan;
aop = orbit.aop;
nu = orbit.nu;

% Central body gravitational constant
switch central_body
    case 'sun'
        mu = mu_sun;
    case 'earth'
        mu = mu_earth;
end

%%--------------------------------------------------------------------------------------
p = sma*(1-ecc ^2);
r_0 = p / (1 + ecc * cos(nu));
%
%%--------------- Coordinates in the perifocal reference system Oxyz -----------------%
%
% position vector coordinates
x = r_0 * cos(nu);
y = r_0 * sin(nu);
%
%
% velocity vector coordinates
Vx_ = -(mu/p)^(1/2) * sin(nu);
Vy_ = (mu/p)^(1/2) * (ecc + cos(nu));
%
%
%%-------------- the geocentric-equatorial reference system OXYZ ---------------------%
%
% position vector components X, Y, and Z
X = (cos(raan) * cos(aop) - sin(raan) * sin(aop) * cos(inc)) * x + (-cos(raan) * sin(aop) - sin(raan) * cos(aop) * cos(inc)) * y;
Y = (sin(raan) * cos(aop) + cos(raan) * sin(aop) * cos(inc)) * x + (-sin(raan) * sin(aop) + cos(raan) * cos(aop) * cos(inc)) * y;
Z = (sin(aop) * sin(inc)) * x + (cos(aop) * sin(inc)) * y;
% velocity vector components X', Y', and Z'
Vx = (cos(raan) * cos(aop) - sin(raan) * sin(aop) * cos(inc)) * Vx_ + (-cos(raan) * sin(aop) - sin(raan) * cos(aop) * cos(inc)) * Vy_;
Vy = (sin(raan) * cos(aop) + cos(raan) * sin(aop) * cos(inc)) * Vx_ + (-sin(raan) * sin(aop) + cos(raan) * cos(aop) * cos(inc)) * Vy_;
Vz = (sin(aop) * sin(inc)) * Vx_ + (cos(aop) * sin(inc)) * Vy_;

% Pack outputs into vectors
position = [X; Y; Z];
velocity = [Vx; Vy; Vz];


