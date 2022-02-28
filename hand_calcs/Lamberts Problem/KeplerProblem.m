function [r2, v2]=KeplerProblem(orbitalelements, deltatime, rvector, vvector)
%mu of sun
Mu=1.3271233e11;
%get the needed orbital elements
a=orbitalelements.sma;
trueanomaly=orbitalelements.nu;
e=orbitalelements.ecc;
%mean motion
meanmotion=(Mu/(a)^3)^(1/2);
%calculate eccentric anomaly
Eccentricanomaly=acos((e+cos(trueanomaly))/(1+e*cos(trueanomaly)));
%initial mean anomaly
meananomaly=-Eccentricanomaly+e*sin(Eccentricanomaly);
%final mean anomaly
Meananomaly2=meananomaly+meanmotion*(deltatime*24*60*60);
%solve keplers equation for final eccentric anomaly
fsolver=@(E) e*sin(E)-E-Meananomaly2;
eccentricanomaly2=fzero(fsolver, Meananomaly2);
%compute the change in eccentric anomaly
deltaE=eccentricanomaly2-Eccentricanomaly;
%compute lagrangian f,g,fdot, gdot
f=1-a/norm(rvector)*(1-cos(deltaE));
g=deltatime-((a)^3/Mu)^(1/2)*(-1*sin(deltaE)+deltaE);
%final position
r2=f*rvector+g*vvector;
fdot=-(Mu*a)^(1/2)*sin(deltaE)/(norm(r2)*norm(rvector));
gdot=1-a/norm(r2)*(1-cos(deltaE));
%final velocity
v2=fdot*rvector+gdot*vvector;

end