function [dataTable] = save_data(sat)
ASTG = sat.Propagator;
%% Identify all target sequences in the MCS
ts_ind = []; % List of indices of target sequences in the MCS
MCS_items = ASTG.MainSequence.Count; % Number of main segments in the MCS
for i = 0:MCS_items-1
    % Check if current segment is a target sequence
    if ASTG.MainSequence.Item(i).Type == "eVASegmentTypeTargetSequence"
        ts = [ts_ind i]; % Store index of target sequence
        
    end
    
end
for jj = 1:i
    ts = ASTG.MainSequence.Item(0); % Gives you the target sequence object
    segment = ts.Segments.Item(0); % Gives you the first item within the target sequence
    % code to check whether "segment" is a maneuver here
    
end
dataTable = 1;
end