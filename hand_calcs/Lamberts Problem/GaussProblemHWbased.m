function [V1,V2]=GaussProblem(rinterceptor,rtargetfinal,deltatknown)
        
    r1=rinterceptor;
    r2=rtargetfinal;
    %short way p iteration
    deltanu=acos(dot(r1,r2)/(norm(r1)*norm(r2)));
    disp(deltanu)
    MU=1.3271233e11;
    k=norm(r1)*norm(r2)*(1-cos(deltanu));
    l=norm(r1)+norm(r2);
    m=norm(r1)*norm(r2)*(1+cos(deltanu));
    peye=k/(l+sqrt(2*m));
    peyeeye=k/(l-sqrt(2*m));
    p0=8926.9
    p1=9373.3


    afunc=@(param) m*k*param/((2*m-l^2)*param^2+2*k*l*param-k^2);
    ffunc=@(param) 1-norm(r2)/param*(1-cos(deltanu));
    fdot=@(param)sqrt(MU/param)*tan(deltanu/2)*((1-cos(deltanu))/param-1/norm(r1)-1/norm(r2));
    cosE=@(param) (1-norm(r1)/(afunc(param))*(1-ffunc(param)));
    sinE=@(param) -norm(r1)*norm(r2)*fdot(param)/sqrt(MU*afunc(param))
    
    %E=@(param) acos(1-norm(r1)/(afunc(param))*(1-ffunc(param)));
    g=@(param) norm(r1)*norm(r2)*sin(deltanu)/sqrt(MU*param);
    gdot=@(param)1-norm(r1)/param*(1-cos(deltanu));
    toffunc=@(param) g(param)+sqrt(afunc(param)^3/MU)*(E(param)-sin(E(param)))
    tofp0=toffunc(p0)
    tofp1=toffunc(p1)
    t1=tofp1
    t0=tofp0
    itnum=0
    dif=6
    while dif>1
        %dt/dp=@(param) -g(param)/(2*param)-3/2*afunc(param)*(deltatknown-g(param)*((k^2+(2*m-l^2)*param^2)/(m*k*p^2))+sqrt(afunc(param)^3/MU)*(2*k*sin(E(param))/(param*(k-l*param))
        pn=p1+(deltatknown-t1)*(p1-p0)/(t1-t0)
        tn=toffunc(pn)
        p0=p1
        p1=pn
        t0=t1
        t1=tn
        dif=deltatknown-t1
        itnum=itnum+1
    end
    %disp(p1)
    %disp(itnum)
    finalf=ffunc(p1);
    finalg=g(p1);
    finalgdot=gdot(p1);
    finalfdot=fdot(p1);
    vx1=(r2(1)-finalf*r1(1))/finalg;
    vy1=(r2(2)-finalf*r1(2))/finalg;
    vz1=(r2(3)-finalf*r1(3))/finalg;
    V1=[vx1;vy1;vz1];
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
