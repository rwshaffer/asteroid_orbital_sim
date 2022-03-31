function [dataTable] = save_data(sat)
ASTG = sat.Propagator;
%% Identify all target sequences in the MCS
ts_ind = []; % List of indices of target sequences in the MCS
MCS_items = ASTG.MainSequence.Count; % Number of main segments in the MCS
deltaV = 0;
for i = 0:MCS_items-1
    item = ASTG.MainSequence.Item(i);
    % Check if current segment is a target sequence
    if item.Type == "eVASegmentTypeTargetSequence"
        countTS = item.Segments.Count;
        %% Check if the segments are maneuvers
        for jj = 0:length(countTS)
            segment = item.Segments.Item(jj); % Gives you the first item within the target sequence
            % code to check whether "segment" is a maneuver here
            isManeuver(i+1,jj+1) = (segment.Type == "eVASegmentTypeManeuver");
            
            if isManeuver(i+1, jj+1)
                deltaV = deltaV + segment.GetResultValue('DeltaV');
            end
        end
    end
    dataTable = deltaV;
end
    disp(['Total deltaV for this Trajectory is: ' num2str(deltaV) 'km/s.'])