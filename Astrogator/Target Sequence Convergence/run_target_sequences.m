function [sat,diverged] = run_target_sequences(sat)


%% Convergence parameters to play with
kmax = 30; % Max number of times to run a given target sequence before giving up
consecutive_divergences_allowed = 1; % Number of RunMCS iterations that any result can diverge for before making changes
div_threshold = 0.5; % Percentage of improvement that a result must see in order to not be considered "diverging"
tolerance_scale_factor = 2; % Factor by which to increase tolerance of a result
max_tolerance_increases = 1; % Number of times it is allowed to increase a result's tolerance before moving on to
                             % more drastic measures
num_IC_changes = 4; % Number of times a TS will have initial conditions randomized if step size increase fails

                             
%% Connection to STK
ASTG = sat.Propagator;

%% Identify all target sequences in the MCS
ts_ind = []; % List of indices of target sequences in the MCS
MCS_items = ASTG.MainSequence.Count; % Number of main segments in the MCS
for i = 0:MCS_items-1
    % Check if current segment is a target sequence
    if ASTG.MainSequence.Item(i).Type == "eVASegmentTypeTargetSequence"
        ts_ind = [ts_ind i]; % Store index of target sequence
    end
end


%% Iterate through all target sequences to get them to converge
for i = ts_ind
    %% Setup STK objects
    ts = ASTG.MainSequence.Item(i);
    dc = ts.Profiles.Item(0); % Assumes each target sequence only has one differential corrector
    ts.Action = 'eVATargetSeqActionRunActiveProfiles'; % Allow current target sequence to be active
    
    %% Track convergence/divergence
    converged = false; % Keep track of convergence
    diverged = false; % Keep track of divergence
        
    %% Keep track of results tolerances and achieved differences
    results_count = dc.Results.Count;
    % Find indices of enabled results
    enabled_results = [];
    for result_ind = 0:results_count-1
        if dc.Results.Item(result_ind).Enable
            enabled_results = [enabled_results, result_ind];
        end
    end
    num_results = length(enabled_results);
    
    tolerances = zeros(1,num_results); % Initialize row vector with result tolerances
    diffs = zeros(1,num_results); % Initialize matrix of results by RunMCS iteration
    percent_changes = 100*ones(1,num_results); % Initialize vector of changes between RunMCS iterations

    % Populate tolerances vector
    for j = 1:num_results
        tolerances(j) = dc.Results.Item(enabled_results(j)).Tolerance;
        
    end
    
    %% Track parameters used to correct target sequence and force convergence
    steps_already_changed = false; % Keep track of whether step sizes have been increased already
    keep_changing_tolerances = true; % Keep track of whether to continue increasing tolerances or give up
    tol_increased = zeros(1,results_count); % Initialize vector used to track whenever a result
                                                    % tolerance is increased due to divergence 
    IC_changed = 0;
    nominal_vals = [];
    IC_scales = [];
    
    disabled_vmag = false;
    vmag_index = nan;
    for j = 1:num_results
        if dc.Results.Item(enabled_results(j)).Name == "V Mag"
            vmag_index = j;
        end
    end
                                                    
    %% Repeatedly run MCS and keep track of results
    k = 0; % Track number of "RunMCS" iterations
    while ~converged && ~diverged && k < kmax
        %% Run MCS, update iteration count, and check for convergence/divergence
        ASTG.RunMCS();
        k = k + 1;
        
        status = dc.Status; % Resulting status of diff. corrector after latest MCS run
        switch status
            case "Converged"
                % First check if a result (like V Mag) was disabled
                if disabled_vmag
                   [dc,disabled_vmag] = enable_vmag(dc); 
                else
                    converged = true;
                    fprintf('Target Sequence %d, Iteration %d: Converged\n\n',i+1,k)
                    break
                end
            case "Encountered an Error"
                diverged = true;
                fprintf('Target Sequence %d, Iteration %d: ERROR\n\n',i+1,k)
                break
        end
        if k >= kmax
            diverged = true;
            fprintf('Maximum number of iterations reached.\n')
            break
        end
        
        %% Keep track of achieved differences on enabled results
        
        for j = 1:num_results
            diffs(k,j) = abs(dc.Results.Item(enabled_results(j)).Difference);
        end

        %% Display some output text to give info on target sequence convergence
        fprintf('Target Sequence %d, Iteration %d:\n',i+1,k)
        fprintf('Percent of Tolerance Achieved:\n')
        disp((diffs(k,:)./tolerances)*100)
        
        %% Determine whether results are improving or if corrective actions are needed
        % NOTE: "Improvement" of a result means that its percent change in difference (which 
        % we want to go to zero) is more negative than -div_threshold, meaning that the results
        % improved by at least div_threshold (defined near the top of the function, around 1%).
        % In other words, a result getting better by 0.1% is not really improving, since we're too
        % impatient to wait for that to converge (and it might never).
        
        if k > 1 % Only check for improvements after the 2nd iteration
            percent_changes(k-1,:) = (diffs(k,:) - diffs(k-1,:))./diffs(k-1,:) * 100; % Vector of changes in each result from one run to the next
            fprintf('Percent Improvements over last iteration:\n')
            disp(-percent_changes(k-1,:))
                
            % Check whether all results improved. If so, apply the changes and skip over next analysis
            if all(percent_changes(k-1,:) < -div_threshold | diffs(k,:)<tolerances)
                ts.ApplyProfileByName("Differential Corrector"); % Assumes diff corr. is named "Differential Corrector"

            % If results did not all improve, more analysis needed:
            else
                % Track whether each result has diverged more than X times in a row
                %   where "X" is the variable div_threshold
                consec_div = track_target_improvement(percent_changes,div_threshold);
                action_needed = (diffs(k,:) >= tolerances) .* ...
                    (consec_div >= consecutive_divergences_allowed);
                %   Note: this vector is 1 for each result if the result is outside its tolerance 
                %   and has diverged too many times in a row (and thus should have some corrective
                %   action taken) and 0 otherwise.
                
                
                %% Handle corrective actions to allow the Target Sequence to converge
                if any(action_needed) 
                    
                    % Check whether V Mag is an enabled result and is converged - if so, disable it
                    if ~isnan(vmag_index) && ~action_needed(vmag_index) && dc.Results.Item(enabled_results(vmag_index)).Enable
                        [ts,dc,disabled_vmag] = remove_vmag(ts,dc);
                    
                    % Check whether step sizes have already been changed
                    elseif ~steps_already_changed
                        % Increase step sizes on duration and epoch stopping conditions within TS
                        [ts,dc] = increase_step_sizes(ts,dc); % Increase relevant step sizes
                        steps_already_changed = true; % Only change step sizes X number of times
                         
                        
                    % Check whether the tolerance increasing function wants to keep trying
                    elseif keep_changing_tolerances
                        % Pack tolerance updating info into a struct to keep code readable
                        tolerance_info.tol_increased = tol_increased;
                        tolerance_info.max_tolerance_increases = max_tolerance_increases;
                        tolerance_info.tolerance_scale_factor = tolerance_scale_factor;
                        % Increase tolerances on results that are struggling to converge
                        [ts,dc,tolerances,tolerance_info,keep_changing_tolerances] = increase_tolerances(ts,dc,tolerances,tolerance_info,action_needed,enabled_results);
                        tol_increased = tolerance_info.tol_increased;
                        
                    % Check whether initial conditions have been changed more than allowed
                    elseif IC_changed < num_IC_changes
                        % Randomize some initial conditions within the target sequence
                        [ts,dc,nominal_vals,IC_scales] = change_init_conditions(ts,dc,IC_changed,nominal_vals,IC_scales);
                        IC_changed = IC_changed + 1; % Update counter on number of times this has been done
   
                        
                    % If all above steps have failed, allow the simulation to end as "diverged"
                    else
                        diverged = true;
                        fprintf('Target Sequence %d DID NOT CONVERGE.\n\n',i+1)
                        break                           
                    end                        
           
                end
              
            end

        end
    
    end
    
    %% Clean up after convergence or divergence
    if converged
        % Apply changes and set target sequence to nominal to speed up future runs
        ts.ApplyProfileByName("Differential Corrector"); % Assumes diff corr. is named "Differential Corrector"
        ts.Action = 'eVATargetSeqActionRunNominalSeq'; % Set converged TS to nominal
    end
    if diverged
        break % We're done with this trajectory, so don't bother with later target sequences
    end    
end






end