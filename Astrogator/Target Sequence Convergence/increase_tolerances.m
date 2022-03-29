function [ts,dc,tolerances,tol_increased,keep_changing_tolerances] = increase_tolerances(ts,dc,tolerances,tol_increased,action_needed,num_tolerance_increases,tolerance_scale_factor)

%fprintf('Action: increasing tolerances\n')
keep_changing_tolerances = true;

results_count = length(action_needed);

for j = 0:results_count-1
    if action_needed(j+1)
        % Check to see if this tolerance has already been increased.
        %   We probably only want to do this once per result
        if tol_increased(j+1) < num_tolerance_increases
            % Update tolerance by a scale factor (defined near top of function) in STK
            dc.Results.Item(j).Tolerance = dc.Results.Item(j).Tolerance*tolerance_scale_factor;
            % Update tolerance in this vector (not sure if it's being used)
            tolerances(j+1) = dc.Results.Item(j).Tolerance;
            % Update vector tracking tolerance increases
            tol_increased(j+1) = tol_increased(j+1) + 1;

            fprintf('Action: increasing tolerance on result %d\n',j+1)
        else
            fprintf('Note: past tolerance increase limit on result %d\n',j+1)
            keep_changing_tolerances = false;
        end

    end

end



end