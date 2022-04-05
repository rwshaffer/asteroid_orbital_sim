function [ts,dc,nominal_vals,IC_scales] = change_init_conditions(ts,dc,IC_changed,nominal_vals,IC_scales)

%% Randomly change some of the variables in the target sequence
% (can be control variables in the diff. corrector, but maybe don't have to be?)
% Look for the following types of variables:
% - Thrust efficiency
% - Impulsive maneuver magnitude
% - Epoch/Duration stopping conditions


num_segments = ts.Segments.Count;

%% Handle the "initial" initial conditions - i.e., the ones hard-coded in
% If this is the first time changing IC's, store IC's as nominal
if IC_changed == 0
    ii = 0; % Number of variables that have been added
    
    for j = 0:num_segments-1
        segment = ts.Segments.Item(j);
        if segment.Type == "eVASegmentTypeManeuver"
            if segment.ManeuverType == "eVAManeuverTypeImpulsive"
                imp_mnvr = segment.Maneuver;
                if imp_mnvr.AttitudeControl.CoordType == "eVASphericalImpDeltaV"
                    % Impulsive burn spherical magnitude will be varied
                    ii = ii + 1;
                    nominal_vals(ii) = imp_mnvr.AttitudeControl.Magnitude;
                    IC_scales(ii) = 0.5; % Allow for varying magnitude by +- 50%
                    
                elseif imp_mnvr.AttitudeControl.CoordType == "eVACartesianImpDeltaV"
                    % Vary impulsive burn Cartesian coordinates
                    ii = ii + 1;
                    nominal_vals(ii) = imp_mnvr.AttitudeControl.X;
                    IC_scales(ii) = 0.5;
                    ii = ii + 1;
                    nominal_vals(ii) = imp_mnvr.AttitudeControl.Y;
                    IC_scales(ii) = 0.5;
                    ii = ii + 1;
                    nominal_vals(ii) = imp_mnvr.AttitudeControl.Z;
                    IC_scales(ii) = 0.5;
                end
            elseif segment.ManeuverType == "eVAManeuverTypeFinite"
                fin_mnvr = segment.Maneuver;
                % Finite maneuver thrust efficiency will be varied
                ii = ii + 1;
                nominal_vals(ii) = fin_mnvr.ThrustEfficiency;
                IC_scales(ii) = 0.5;
                % Duration stopping condition will be varied
                stopping_cond = fin_mnvr.Propagator.StoppingConditions.Item(0);
                if stopping_cond.Name == "Duration"
                    ii = ii + 1;
                    nominal_vals(ii) = stopping_cond.Properties.Trip;
                    IC_scales(ii) = 0.1; % Vary duration by up +- 10%
                % Epoch stopping condition will be varied
                elseif stopping_cond.Name == "Epoch"
                    ii = ii + 1;
                    % Store the value as a Julian Date (for now)
                    date = stopping_cond.Properties.Trip;
                    dateFormatString = 'd MMM yyyy H:mm:ss.SSS';
                    nominal_vals(ii) = juliandate(datetime(date, 'InputFormat', dateFormatString));
                    IC_scales(ii) = 0.00002;
                end
                
            end
        elseif segment.Type == "eVASegmentTypePropagate"
            % Duration stopping condition will be varied
            stopping_cond = segment.StoppingConditions.Item(0);
            if stopping_cond.Name == "Duration"
                ii = ii + 1;
                nominal_vals(ii) = stopping_cond.Properties.Trip;
                IC_scales(ii) = 0.1; % Vary duration by up +- 10%
            % Epoch stopping condition will be varied
            elseif stopping_cond.Name == "Epoch"
                ii = ii + 1;
                % Store the value as a Julian Date (for now)
                date = stopping_cond.Properties.Trip;
                dateFormatString = 'd MMM yyyy H:mm:ss.SSS';
                nominal_vals(ii) = juliandate(datetime(date, 'InputFormat', dateFormatString));
                IC_scales(ii) = 0.00002;
            end
        end
    end
else
    ii = length(nominal_vals);
end


%% Randomize a random set of initial conditions
vars_to_vary = randi([0,1],1,ii); % Vector of 1s and 0s, where a 1 means that the variable 
% gets changed from nominal, and a 0 means it remains nominal

for j = 1:ii
    if vars_to_vary(j)
        new_vals(j) = nominal_vals(j)*(1 + IC_scales(j)*(-1+2*rand));
    else
        new_vals(j) = nominal_vals(j);
    end
end


%% Iterate through and randomize conditions
% Since the for loop will go through the segments in the same order and
% look for the same things, the indices ii still match up with the correct
% nominal values/scale factors for each variable
ii = 0; % Reset counter

for j = 0:num_segments-1
    segment = ts.Segments.Item(j);
    if segment.Type == "eVASegmentTypeManeuver"
        if segment.ManeuverType == "eVAManeuverTypeImpulsive"
            imp_mnvr = segment.Maneuver;
            if imp_mnvr.AttitudeControl.CoordType == "eVASphericalImpDeltaV"
                % Impulsive burn spherical magnitude will be varied
                ii = ii + 1;
                imp_mnvr.AttitudeControl.Magnitude = new_vals(ii);

            elseif imp_mnvr.AttitudeControl.CoordType == "eVACartesianImpDeltaV"
                % Vary impulsive burn Cartesian coordinates
                ii = ii + 1;
                imp_mnvr.AttitudeControl.X = new_vals(ii);
                ii = ii + 1;
                imp_mnvr.AttitudeControl.Y = new_vals(ii);
                ii = ii + 1;
                imp_mnvr.AttitudeControl.Z = new_vals(ii);
            end
        elseif segment.ManeuverType == "eVAManeuverTypeFinite"
            fin_mnvr = segment.Maneuver;
            % Finite maneuver thrust efficiency will be varied
            ii = ii + 1;
            fin_mnvr.ThrustEfficiency = new_vals(ii);
            % Duration stopping condition will be varied
            stopping_cond = fin_mnvr.Propagator.StoppingConditions.Item(0);
            if stopping_cond.Name == "Duration"
                ii = ii + 1;
                stopping_cond.Properties.Trip = new_vals(ii);
            % Epoch stopping condition will be varied
            elseif stopping_cond.Name == "Epoch"
                ii = ii + 1;
                % Store the value as a Julian Date (for now)
                % STILL NEED TO FIGURE THIS PART OUT!!!
                %stopping_cond.Properties.Trip = nominal_vals(ii);
                
            end

        end
    elseif segment.Type == "eVASegmentTypePropagate"
        % Duration stopping condition will be varied
        stopping_cond = segment.StoppingConditions.Item(0);
        if stopping_cond.Name == "Duration"
            ii = ii + 1;
            stopping_cond.Properties.Trip = new_vals(ii);
        % Epoch stopping condition will be varied
        elseif stopping_cond.Name == "Epoch"
            ii = ii + 1;
            % Store the value as a Julian Date (for now)
            % STILL NEED TO FIGURE THIS PART OUT!!!
            %stopping_cond.Properties.Trip = nominal_vals(ii);
        end
    end
end


% Reset target sequence to last time it had all results improve
ts.ResetProfileByName("Differential Corrector");



fprintf('Action: randomly changing some initial conditions\n\n')






end