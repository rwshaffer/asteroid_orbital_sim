function [V1,V2]=GaussProblemtextbook(rinterceptor,rtargetfinal,deltatknown)
try
    %convert from days to seconds
    format long
    deltatknown;
    deltatknown=deltatknown*24*60*60;
    r1=reshape(rinterceptor,3,1);
    r2=reshape(rtargetfinal, 3, 1);
    %short way p iteration
    deltanu=acos(dot(r1,r2)/(norm(r1)*norm(r2)));
    %disp(deltanu)
    MU=1.3271233e11;
    k=norm(r1)*norm(r2)*(1-cos(deltanu));
    l=norm(r1)+norm(r2);
    m=norm(r1)*norm(r2)*(1+cos(deltanu));
    peye=k/(l+sqrt(2*m));
    peyeeye=k/(l-sqrt(2*m));
    p1=(peye+peyeeye)/2;


    afunc=@(param) m*k*param/((2*m-l^2)*param^2+2*k*l*param-k^2);
    ffunc=@(param) 1-norm(r2)/param*(1-cos(deltanu));
    fdot=@(param)sqrt(MU/param)*tan(deltanu/2)*((1-cos(deltanu))/param-1/norm(r1)-1/norm(r2));
    cosE=@(param) (1-norm(r1)/(afunc(param))*(1-ffunc(param)));
    sinE=@(param) -norm(r1)*norm(r2)*fdot(param)/sqrt(MU*afunc(param));
    E=@(param) atan2(sinE(param), cosE(param));
    %E=@(param) acos(1-norm(r1)/(afunc(param))*(1-ffunc(param)));
    g=@(param) norm(r1)*norm(r2)*sin(deltanu)/sqrt(MU*param);
    gdot=@(param)1-norm(r1)/param*(1-cos(deltanu));
    toffunc=@(param) g(param)+sqrt(afunc(param)^3/MU)*(E(param)-sin(E(param)));
    tn=toffunc(peyeeye);
    itnum=0;
    tol=10^-2;
    dt_dp=@(param) -g(param)/(2*param)-3/2*afunc(param)*(deltatknown-g(param))*((k^2+(2*m-l^2)*param^2)/(m*k*param^2))+sqrt(afunc(param)^3/MU)*(2*k*sin(E(param))/(param*(k-l*param)));
        
    while abs(tn-deltatknown)>tol
        derivative=dt_dp(p1);
        pn=p1+(deltatknown-tn)/derivative;
        tn=toffunc(pn);
        p1=pn;
        itnum=itnum+1;
        if itnum>10000
            error("Iterations exceeded, will not converge.")
        end
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
    pause
catch
    %this makes this an unacheivable delta v
    disp("Gauss Failed")
    V1=[1000000000; 0; 0];
    V2=[0;0;0];
end

end
