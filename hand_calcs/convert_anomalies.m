function [nu,E,M] = convert_anomalies(input,ecc,in_type)
% Convert between true, eccentric, and mean anomalies.
% ALL ANGLES IN RADIANS!
%
% INPUTS
% input - input anomaly (radians)
% ecc - orbit eccentricity
% in_type - "true", "ecc", or "mean" - starting anomaly type
%
% OUTPUTS
% nu - true anomaly (radians)
% E - eccentric anomaly (radians)
% M - mean anomaly (radians)

switch in_type
%% Convert from true anomaly to eccentric or mean anomaly
    case "true"
        nu = input;
        E = acos((ecc + cos(nu))/(1 + ecc*cos(nu)));
        if nu > pi
            E = 2*pi - E;
        end
        M = E - ecc * sin(E); % Quad check or no?
%% Convert from eccentric anomaly to true or mean anomaly        
    case "ecc"
        E = input;
        nu = acos((cos(E) - ecc)/(1 - ecc*cos(E)));
        if E > pi
            nu = 2*pi - nu;
        end
        M = E - ecc * sin(E);
%% Convert from mean anomaly to true or eccentric anomaly
    case "mean"
        M = input;
        E = CalcEA(M,ecc); % Calls separate function for Newton-Raphson solution
        nu = atan2((sin(E)*(1-ecc^2)^.5),(cos(E)-ecc)); % Double check this. Quad check?
        if E > pi
            nu = 2*pi + nu;
        end
        
    otherwise
        disp('Invalid anomaly type.')
end


end