function [dataTable] = extract_data(sat,mining_duration,outbound_traj,Isp)
%% Handle optional Isp input
if nargin < 4
    Isp = 1710;
end


ASTG = sat.Propagator;
%% Outbound Delta-V
MCS_items = ASTG.MainSequence.Count; % Number of main segments in the MCS

outbound_fin_deltaV = 0; % Finite delta-V on outbound trajectory
outbound_imp_deltaV = 0; % Impulsive delta-V on outbound trajectory
max_thrust = 0; % Maximum thrust required over trajectory
for i = 0:MCS_items-3 % Hack: save the last target sequence (return traj. for later)
    item = ASTG.MainSequence.Item(i);
    % Check if current segment is a target sequence
    if item.Type == "eVASegmentTypeTargetSequence"
        countTS = item.Segments.Count;
        %% Check if the segments are maneuvers
        for jj = 0:countTS-1
            segment = item.Segments.Item(jj); % Gives you the first item within the target sequence
            isManeuver = (segment.Type == "eVASegmentTypeManeuver");
            
            if isManeuver
                if segment.ManeuverType == "eVAManeuverTypeFinite"
                    outbound_fin_deltaV = outbound_fin_deltaV + segment.GetResultValue('DeltaV');
                    % Check if this finite burn uses the max. thrust
                    segment_thrust = segment.Maneuver.ThrustEfficiency;
                    if segment_thrust > max_thrust
                        max_thrust = segment_thrust;
                    end
                elseif segment.ManeuverType == "eVAManeuverTypeImpulsive"
                    outbound_imp_deltaV = outbound_imp_deltaV + segment.GetResultValue('DeltaV');
                end
            end
        end
    % Do the same thing for any backwards sequence, but first find the TS inside
    elseif item.Type == "eVASegmentTypeBackwardSequence"
        bs_ind = i; % MCS index of backwards sequence
        item = item.Segments.Item(0);
        countTS = item.Segments.Count;
        %% Check if the segments are maneuvers
        for jj = 0:countTS-1
            segment = item.Segments.Item(jj); % Gives you the first item within the target sequence
            isManeuver = (segment.Type == "eVASegmentTypeManeuver");
            
            if isManeuver
                if segment.ManeuverType == "eVAManeuverTypeFinite"
                    % Note: this uses some janky stuff to get deltaV, since it seems that STK gives backwards finite maneuver delta V results as 0
                    % Get S/C mass at beginning of burn (end of segment)
                    seg_init_mass = segment.FinalState.DryMass + segment.FinalState.FuelMass;
                    % Get S/C mass at end of burn (beginning of segment)
                    seg_fin_mass = segment.InitialState.DryMass + segment.InitialState.FuelMass;
                    % Apply the rocket equation to get segment delta V in km/s
                    seg_deltaV = Isp*9.81*log(seg_init_mass/seg_fin_mass) * 1e-3;
                    % Add this delta V to the running sum
                    outbound_fin_deltaV = outbound_fin_deltaV + seg_deltaV;
                    % Check if this finite burn uses the max. thrust
                    segment_thrust = segment.Maneuver.ThrustEfficiency;
                    if segment_thrust > max_thrust
                        max_thrust = segment_thrust;
                    end
                elseif segment.ManeuverType == "eVAManeuverTypeImpulsive"
                    outbound_imp_deltaV = outbound_imp_deltaV + segment.GetResultValue('DeltaV');
                end
            end
        end
    end
    
end

%% Return Delta-V
% return trajectory is always the last target sequence
return_deltaV = 0; % Note that there are no impulsive burns on return traj.
i = i + 1;
item = ASTG.MainSequence.Item(i);
if item.Type == "eVASegmentTypeTargetSequence"
    countTS = item.Segments.Count;
    %% Check if the segments are maneuvers
    for jj = 0:countTS-1
        segment = item.Segments.Item(jj); % Gives you the first item within the target sequence
        isManeuver = (segment.Type == "eVASegmentTypeManeuver");

        if isManeuver
            return_deltaV = return_deltaV + segment.GetResultValue('DeltaV');
            % Check if this finite burn uses the max. thrust
            if segment.ManeuverType == "eVAManeuverTypeFinite"
                segment_thrust = segment.Maneuver.ThrustEfficiency;
                if segment_thrust > max_thrust
                    max_thrust = segment_thrust;
                end
            end
        end
    end
end

max_thrust = max_thrust * 100; % Thrust efficiency in percent

%% C3 Energy Required
% Note: this will play into delta-V required, but we'll handle it during data processing
% The current approach only works for direct outbound trajectories.
switch outbound_traj
    case "Direct"
        initial_state = ASTG.MainSequence.Item(0).Segments.Item(0);
        C3_Energy = initial_state.Element.C3Energy;
        C3_Energy = max(C3_Energy,0); % If C3 comes out negative, set it to 0
    case "EGA"
        back_seq = ASTG.MainSequence.Item(bs_ind);
        back_seq_ts = back_seq.Segments.Item(0);
        ts_count = back_seq_ts.Segments.Count;
        final_segment = back_seq_ts.Segments.Item(ts_count-2);
        C3_Energy = final_segment.GetResultValue('C3_Energy');
        C3_Energy = max(C3_Energy,0); % If C3 comes out negative, set it to 0
end
%% Outbound Time of Flight
switch outbound_traj
    case "Direct"
        first_ts = ASTG.MainSequence.Item(0);
        first_segment = first_ts.Segments.Item(0);
        t0_outbound = first_segment.InitialState.Epoch;

        num_ts = ASTG.MainSequence.Count;
        final_outbound_ts = ASTG.MainSequence.Item(num_ts - 3);
        num_segments = final_outbound_ts.Segments.Count;
        final_outbound_segment = final_outbound_ts.Segments.Item(num_segments-2);
        tf_outbound = final_outbound_segment.FinalState.Epoch;
    case "EGA"
        back_seq = ASTG.MainSequence.Item(bs_ind);
        back_seq_ts = back_seq.Segments.Item(0);
        ts_count = back_seq_ts.Segments.Count;
        final_segment = back_seq_ts.Segments.Item(ts_count-2);
        t0_outbound = final_segment.FinalState.Epoch;
        
        num_ts = ASTG.MainSequence.Count;
        final_outbound_ts = ASTG.MainSequence.Item(num_ts - 4);
        num_segments = final_outbound_ts.Segments.Count;
        final_outbound_segment = final_outbound_ts.Segments.Item(num_segments-2);
        tf_outbound = final_outbound_segment.FinalState.Epoch;
end

%process t0 and tf into MATLAB datetime objects
dateFormatString = 'd MMM yyyy H:mm:ss.SSS';
t0_outbound = datetime(t0_outbound, 'InputFormat', dateFormatString);
tf_outbound = datetime(tf_outbound, 'InputFormat', dateFormatString);

%calculate elapsed time
ToF_outbound = tf_outbound - t0_outbound;
ToF_outbound.Format = 'd';

%format elapsed time into string for output
ToF_String_outbound = char(ToF_outbound, 'd');

%% Return Time of Flight
last_ts = ASTG.MainSequence.Item(MCS_items-2);
first_segment = last_ts.Segments.Item(0);
t0_return = first_segment.InitialState.Epoch;

num_segments = last_ts.Segments.Count;
final_return_segment = last_ts.Segments.Item(num_segments-2);
tf_return = final_return_segment.FinalState.Epoch;

%process t0 and tf into MATLAB datetime objects
dateFormatString = 'd MMM yyyy H:mm:ss.SSS';
t0_return = datetime(t0_return, 'InputFormat', dateFormatString);
tf_return = datetime(tf_return, 'InputFormat', dateFormatString);

%calculate elapsed time
ToF_return = tf_return - t0_return - days(mining_duration);
ToF_return.Format = 'd';

%format elapsed time into string for output
ToF_String_return = char(ToF_return, 'd');

%% Format data for output and display
dataTable = [outbound_fin_deltaV; outbound_imp_deltaV; return_deltaV; max_thrust; C3_Energy; string(t0_outbound,dateFormatString); days(ToF_outbound); days(ToF_return)];
disp(['Total Finite Delta-V for outbound trajectory is: ' num2str(outbound_fin_deltaV) ' km/s.'])
disp(['Total Impulsive Delta-V for outbound trajectory is: ' num2str(outbound_imp_deltaV) ' km/s.'])
disp(['Total Delta-V for return trajectory is: ' num2str(return_deltaV) ' km/s.'])
disp(['Maximum thrust over entire trajectory is: ' num2str(max_thrust) ' %.'])
disp(['C3 Energy required at Earth escape: ' num2str(C3_Energy) ' km^2/s^2.'])
disp(['Date of Earth escape: ' string(t0_outbound,dateFormatString)])
disp(['Time of flight for outbound trajectory: ' ToF_String_outbound '.'])
disp(['Time of flight for return trajectory: ' ToF_String_return '.'])

end