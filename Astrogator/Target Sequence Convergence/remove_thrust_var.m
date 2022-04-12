function [ts,dc] = remove_thrust_var(ts,dc)

%% Find and disable any Thrust Efficiency controls in a target sequence if they are below a threshold
disable_threshold = 0.1;

controls_count = dc.ControlParameters.Count;

for j = 0:controls_count-1
    if dc.ControlParameters.Item(j).Name == "FiniteMnvr.Thrusting.ThrustEfficiency" && dc.ControlParameters.Item(j).Enable
        if dc.ControlParameters.Item(j).FinalValue < disable_threshold
            dc.ControlParameters.Item(j).Enable = 0;
            fprintf("Disabling Control %d, Thrust Efficiency\n\n",j)

        end
        
    end
end


end