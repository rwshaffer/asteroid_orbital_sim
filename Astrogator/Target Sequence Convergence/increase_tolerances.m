function [ts,dc,keep_changing_tolerances] = increase_tolerances(ts,dc,tol_increased)

fprintf('Action: increasing tolerances\n')
keep_changing_tolerances = true;

% for j = 0:results_count-1
%     if action_needed(j+1)
%         % Check to see if this tolerance has already been increased.
%         %   We probably only want to do this once per result
%         if tol_increased(j+1) < num_tolerance_increases
%             % Update tolerance by a scale factor (defined near top of function) in STK
%             dc.Results.Item(j).Tolerance = dc.Results.Item(j).Tolerance*tolerance_scale_factor;
%             % Update tolerance in this vector (not sure if it's being used)
%             tolerances(j+1) = dc.Results.Item(j).Tolerance;
%             % Update vector tracking tolerance increases
%             tol_increased(j+1) = tol_increased(j+1) + 1;
% 
%             fprintf('Increasing tolerance on result %d\n',j+1)
%         else
%             % Apply next corrective action: change step sizes
%         end
% 
%         % Reset target sequence to last time it had its changes applied
%         ts.ResetProfileByName("Differential Corrector");
% 
%     end
% 
% end
end