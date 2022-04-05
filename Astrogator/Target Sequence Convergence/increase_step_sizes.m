function [ts,dc] = increase_step_sizes(ts,dc)
% Look for Control Parameters in the Target Sequence that can have their step size increased
% (based solely on previous experience - this action can improve convergence for some reason).
% The specific parameters that can be scaled are specified in the code below, but are by default
% Epoch and Duration stopping conditions for Propagate/Finite Maneuver segments.
% 
% INPUTS:
% ts: Target Sequence object
% dc: Differential Corrector object within Target Sequence
%
% OUTPUTS:
% ts: Updated Target Sequence object
% dc: Updated Differential Corrector object within Target Sequence


fprintf('Action: increasing step sizes for any Duration/Epoch controls\n\n')

%% List of control parameters to increase step sizes for, and corresponding scaling factors
params = ["Epoch.TripValue","Duration.TripValue"];
scale_factors = [5,5];

%% Find Epoch and Duration stopping conditions among the control parameters
controls_count = dc.ControlParameters.Count; % All parameters used in target sequence

for i = 0:controls_count-1
    control = dc.ControlParameters.Item(i); % i'th control parameter
    control_name = control.Name; % Name of control parameter (including long intro string)
    % Note: e.g. control_name = 'FiniteMnvr.Propagator.StoppingConditions.Epoch' or something
    name_length = length(control_name); % Get length of full name
    
    % Iterate through parameters of interest to see if this control parameter is one of them
    for j = 1:length(params) 
        
        param_word_length = length(char(params(j))); % Length of word at end of control_name to check for
        
        % Check whether the ending string of the control_name is the parameter name
        if control_name(name_length-(param_word_length-1):end) == params(j)
            % Scale the control's max step based on the specified scaling factor
            control.MaxStep = control.MaxStep * scale_factors(j);
        end
    end
end

end