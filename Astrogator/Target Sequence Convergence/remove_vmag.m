function [ts,dc,disabled_vmag] = remove_vmag(ts,dc)

%% Find and disable any V Mag or C3 results in a target sequence

results_count = dc.Results.Count;

for j = 0:results_count-1
    if dc.Results.Item(j).Name == "V Mag" && dc.Results.Item(j).Enable
        dc.Results.Item(j).Enable = 0;
        
        disabled_vmag = true;
        fprintf("Disabling Result %d, V Mag\n\n",j)
        
    elseif dc.Results.Item(j).Name == "C3 Energy" && dc.Results.Item(j).Enable
        dc.Results.Item(j).Enable = 0;
        
        disabled_vmag = true;
        fprintf("Disabling Result %d, C3 Energy\n\n",j)
        
        
    end
end

% Reset target sequence to last time it had all results improve
ts.ResetProfileByName("Differential Corrector");


end