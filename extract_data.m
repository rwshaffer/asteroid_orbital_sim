function [dataTable] = extract_data(sat)
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
        for jj = 0:countTS-1
            segment = item.Segments.Item(jj); % Gives you the first item within the target sequence
            % code to check whether "segment" is a maneuver here
            isManeuver = (segment.Type == "eVASegmentTypeManeuver");
            
            if isManeuver
                deltaV = deltaV + segment.GetResultValue('DeltaV');
            end
        end
    end
    
end
ts = ASTG.MainSequence.Item(0);
first_segment = ts.Segments.Item(0);
t0 = first_segment.InitialState.Epoch;

num_ts = ASTG.MainSequence.Count;
final_ts = ASTG.MainSequence.Item(num_ts - 2);
num_segments = final_ts.Segments.Count;
final_segment = final_ts.Segments.Item(num_segments-2);
tf = final_segment.FinalState.Epoch;

%process t0 and tf into MATLAB datetime objects
dateFormatString = 'd MMM yyyy H:mm:ss.SSS';
t0 = datetime(t0, 'InputFormat', dateFormatString);
tf = datetime(tf, 'InputFormat', dateFormatString);

%calculate elapsed time
ToF = tf - t0;
ToF.Format = 'd';

%format elapsed time into string for output
ToF_String = char(ToF, 'd');

%format data for output
dataTable = [deltaV; days(ToF)];
disp(['Total Delta-V for this trajectory is: ' num2str(deltaV) 'km/s.'])
disp(['The time of flight for this trajectory is: ' ToF_String '.'])
end