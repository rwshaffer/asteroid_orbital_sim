function [V1,V2]=GaussProblem(rinterceptor,rtargetfinal,deltatknown)
        
    r1=rinterceptor;
    r2=rtargetfinal;
    %short way p iteration
    deltanu=dot(r1,r2)/(norm(r1)*norm(r2));
    MU=1.3271233e11;
    k=norm(r1)*norm(r2)*(1-cos(deltanu));
    l=norm(r1)+norm(r2);
    m=norm(r1)*norm(r2)*(1+cos(deltanu));
    peye=k/(l+sqrt(2*m));
    peyeeye=k/(l-sqrt(2*m));
    p0=peye;
    p1=peyeeye;


    afunc=@(param) m*k*param/((2*m-l^2)*param^2+2*k*l*param-k^2);
    ffunc=@(param) 1-norm(r2)/param*(1-cos(deltanu));
    g=@(param) norm(r1)*norm(r2)*sin(deltanu)/sqrt(MU*param);
    gdot=@(param)1-norm(r1)/param*(1-cos(deltanu));
    fdot=@(param)sqrt(MU/param)*tan(deltanu/2)*((1-cos(deltanu))/param-1/norm(r1)-1/norm(r2));
    sinE=@(param) -norm(r1)*norm(r2)*fdot(param)/sqrt(MU*afunc(param));
    cosE=1-(1-ffunc(param))*norm(r1)/a;
    
    E=@(param) atan2(sinE,cosE);
    toffunc=@(param) g(param)+sqrt(afunc(param)^3/MU)*(E(param)-sin(E(param)));
    F0=toffunc(p0)-deltatknown;
    F1=toffunc(p1)-deltatknown;
    itnum=0;
    tol=10^-5;
    %while loop to find orbit based on given time
    while abs(F1)>tol
        deltap=-(p0-p1)*F1/(F0-F1);
        pn=p1+deltap/(1+deltap^2);
        Fn=toffunc(pn)-deltatknown;
        p0=p1;
        p1=pn;
        F0=F1;
        F1=Fn;
        itnum=itnum+1;
    end
    %recompute the lagrangian parameters
    finalf=ffunc(p1);
    finalg=g(p1);
    finalgdot=gdot(p1);
    finalfdot=fdot(p1);
    %find the initial velocity needed for the interceptor
    V1=(r2-finalf.*r1)./finalg;
    %find the final velocity when the interceptor arrives at the target
    V2=finalfdot*r1+finalgdot*V1;
    norm(V1);
    %{
    catch
    %this makes this an unacheivable delta v
    disp("Gauss Failed")
    V1=[1000000000; 0; 0];
    V2=[0;0;0];
    end
%}
end
