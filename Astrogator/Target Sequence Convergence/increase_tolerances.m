function [ts,dc,tolerances,tolerance_info,keep_changing_tolerances] = increase_tolerances(ts,dc,tolerances,tolerance_info,action_needed)

%% Unpack inputs from struct
tol_increased = tolerance_info.tol_increased;
max_tolerance_increases = tolerance_info.max_tolerance_increases;
tolerance_scale_factor = tolerance_info.tolerance_scale_factor;


%% Setup
keep_changing_tolerances = true; % Initialize flag to continue this corrective action vs. move onto the next one
results_count = length(action_needed); % Number of results in the TS to look at

%% Iterate through results
for j = 0:results_count-1
    
    % Check whether the result has been marked as diverging (thus needing action) 
    if action_needed(j+1)
        
        % Check to see if this tolerance has already been increased more than the 
        %   specified max number of tolerance increases
        if tol_increased(j+1) < max_tolerance_increases
            
            % Update tolerance by a scale factor in STK
            dc.Results.Item(j).Tolerance = dc.Results.Item(j).Tolerance*tolerance_scale_factor;
            % Update tolerance in this vector
            tolerances(j+1) = dc.Results.Item(j).Tolerance;
            % Update vector tracking increases in tolerances
            tol_increased(j+1) = tol_increased(j+1) + 1;

            % Display action
            fprintf('Action: increasing tolerance on result %d\n',j+1)
        else
            % This tolerance has been increased too many times. Move on from tolerance increase method
            fprintf('Note: past tolerance increase limit on result %d\n',j+1)
            keep_changing_tolerances = false; % Stop increasing the tolerances, it's not working
        end

    end

end

%% Pack outputs into struct
tolerance_info.tol_increased = tol_increased;


end