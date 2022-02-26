function [deltavdate]=Interceptproblem(MinMaxStartDate, MinMaxtime, interceptororbit, targetorbit)
%{
To do: check anomaly conversions, check gauss problem, retry everything,
debug
%}


% Description of how the code will operate

% will intercept the spacecraft with 1989 ML on a direct trajectory using 
% Impulsive maneuvers. 


% Inputs: A time vector: 
% MinMaxStartDate: Julian Date for beginning of launch window and end of
% launch window
%MinMaxtime: [minimumtime, maximumtime] %given in days
% a structure with inital satellite heliocentric orbit parameters including epoch
% a structure with initial 1989 ML orbit parameters including epoch


% Outputs: The initial delta-V (trajectory starting delta V),final delta-V
% (target incercepting delta-V), and time of flight will all be returned

if nargin< 4
    %define targetMLorbit
    %J2000 coordinates
    targetorbit.epoch=2459600.5; %julian date
    targetorbit.sma=1.9033*10^08; %km
    targetorbit.ecc=0.1363639652278715;
    targetorbit.inc=4.378033480747216;
    targetorbit.raan=104.33806293*pi/180; %rad
    targetorbit.aop=183.360697515989*pi/180; %rad
    targetorbit.nu=3.9826; %rad 
    if nargin<3
        %define satellite orbit
        %will use Earth orbit as approximations (realtively accurate
        %J2000 coordinates
        interceptororbit.epoch=2451545.0; %julian date
        interceptororbit.sma=149.598*10^6; %km
        interceptororbit.ecc=0.01671022;
        interceptororbit.inc=0.00005;
        interceptororbit.raan=3.4542; %rads
        interceptororbit.aop=4.6679; %rads
        interceptororbit.nu=6.1948; %rads
        if nargin<2
            %minmaxtimevector given in days
            MinMaxtime=[50, 1826.25];
            if nargin<1
                %launch dates in julian caldar
                %will use May 2023- Jan 2028 as dates
                MinMaxStartDate=[2460065.5,2461771.5];
            end
        end
    end
end

%% 




date=linspace(MinMaxStartDate(1),MinMaxStartDate(2), 60);
deltavdate=zeros(7,length(date));
for n=1:length(date)
    t=linspace(MinMaxtime(1),MinMaxtime(2), 60);
    deltav=zeros(length(t),1);
    %get epoch in cartestian coordinates
    [rtargetinitial,vtargetinitial]=Kepler2Carts(targetorbit, "sun");
    [rinteceptorepoch,vinterceptorepoch]=Kepler2Carts(interceptororbit, "sun");
    %update satellite/Earth position to launch date
    [rinterceptorinitial,vinterceptorinitial]=KeplerProblem(interceptororbit, date(n)-interceptororbit.epoch,rinteceptorepoch,vinterceptorepoch);
    for i=1:length(t)
        %update position of 1989 ML to after flight, needs date to update from epoch
        [rtargetfinal,vtargetfinal]=KeplerProblem(targetorbit, date(n)-targetorbit.epoch+t(i),rtargetinitial,vtargetinitial);
        %Note: gauss problem will fail for many times
        [v1,v2]=GaussProblemtextbook(rinterceptorinitial,rtargetfinal,t(i));
        deltavinitial=abs(norm(v1-vinterceptorinitial));
        deltavfinal=abs(norm(vtargetfinal-v2));
        deltav(i)=deltavinitial+deltavfinal;
    end
    %find minimum deltav's
    [~,index]=min(deltav);
    %this will be structured to give the 1st row date, 2nd row minimum delta v
    %3rd row time associated, 4th and 5th row one interval below that,
    %6th and 7th row one interval above the minimum. This will ensure the
    %true minimum is within the intervale
    disp(index);
    disp(index(1));
    disp(deltav);
    if index(1)==length(deltav)
        deltavdate(:,n)=[date(n); deltav(index(1)); t(index(1));deltav(index(1)-1); t(index(1)-1);NaN; NaN];
    elseif index(1)==1
        deltavdate(:,n)=[date(n); deltav(index(1)); t(index(1));NaN; NaN;deltav(index(1)+1); t(index(1)+1)];
    else
        deltavdate(:,n)=[date(n); deltav(index(1)); t(index(1));deltav(index(1)-1); t(index(1)-1);deltav(index(1)+1); t(index(1)+1)];
    end
%need to analyze the best launch dates to find what works best for
%our system- then rerun this program using the dates that are selected
%may need to refine the times further to optimize
%may also need to eliminate some options for launch dates based on schedule
end
end
    
%{
probably dont want to use this
    finestep=71.0500*2/100; %gives 101 iterations
    %iterate through the coarse range of optimal times
    tcoarsefine=[tcoarsecoarse(index-1):finestep:tcoarsecoase(index+1)];
    deltavcoarsefine=zeros(length(tcoarsefine),1);
    for i=1:length(tcoarsefine)
        [rtargetfinal,vtargetfinal]=KeplerProblem(rtargetinitial,vtargetinitial,tcoarsefine(i));
        %Note: gauss problem will fail for many times
        [v1,v2]=GaussProblem(rinteceptor,rtargetfinal,tcoarsefine(i));
        deltavinitial=abs(v1-vinterceptor);
        deltavfinal=abs(vtargetfinal-v2);
        deltavcoarsefine(i)=deltavinitial+deltavfinal;
    end
    [~,index]=min(deltavcoarsefine);
    
    deltavcoarsedate(:,n)=[date(n); deltavcoarsefine(index); tcoarsefine(index);deltavcoarsefine(index-1); tcoarsefine(index-1);deltavcoarsefine(index+1); tcoarsefine(index+1)];
    
end

%}

%{
 This is a description of what the code should do
Assuming a restricted two-body problem (neglecting the mass of the satellites
 and no other forces at work) the problem can be solved with search combined 
with a solver for the Gauss' problem.

The gauss problem is the following: Given two positions and the time between
 the positions, find the velocites at both positions, which in other words 
finds the orbit. We will denote this as (v1,v2)=GaussProblem(r1,r2,t).

So how can we use the solver to intercept the second satellite? This is done
by assuming an intercept time (t) (we will iterate to find the best one later).
So what we do now is that we calculates the position of the second satellite at 
the intercept time, which I assume you already can do (using numerical integration 
or solver for the Kepler problem). I will denote this as (r′,v′)=KeplerProblem(r,v,t).

Now that we have two positions (the current position of the first satellite,
the position of the second satellite at the intercept time), we can use the 
solver for the Gauss problem to find what the velocity at the first position 
should have been, if the first satellite was on the intercept orbit.

So our maneuver is simply then: Δv=v1−vsatellite1.

The algorithm is the following:

For each t∈[min,max] (take some range)
(r′,_)=KeplerProblem(rsatellite2,vsatellite2,t)
(v1,v2)=GaussProblem(rsatellite1,r′,t)
Δv=v1−vsatellite1
Select the maneuver that best matches your parameters (minimum intercept time 
given delta velocity requirements).
Some notes:

For a given problem instance, a solver for the gauss problem may fail.
The times are relative, not absolute.
The maneuvers are assumed to be executed at the current time. Calculating 
the intercept at a later time may yield a better intercept in terms of time and deltaV.
Gauss Problem

I'm not 100% sure which is the correct name for the problem, as in 
Fundamentals of astrodynamics by Bate, Roger R., Donald D. Mueller, and Jerry E. White 
denotes this as the Gauss problem, others as Lambert's problem.

For a method for solving the problem, you should either look in the mentioned book,
 which contains several methods, or look at the following links: 
http://aerospacengineering.net/?p=1614, 
http://www.dept.aoe.vt.edu/~cdhall/courses/aoe4134/Apiteration.pdf, 
https://en.m.wikipedia.org/wiki/Lambert%27s_problem.
%}