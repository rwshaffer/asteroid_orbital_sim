function [dc,disabled_vmag] = enable_vmag(dc)

%% Find and disable any V Mag results in a target sequence

results_count = dc.Results.Count;

for j = 0:results_count-1
    if dc.Results.Item(j).Name == "V Mag" && ~dc.Results.Item(j).Enable
        dc.Results.Item(j).Enable = 1;
        
        disabled_vmag = false;
        fprintf("Enabling Result %d, V Mag\n\n",j)
        
    elseif dc.Results.Item(j).Name == "C3 Energy" && ~dc.Results.Item(j).Enable
        dc.Results.Item(j).Enable = 1;
        
        disabled_vmag = false;
        fprintf("Enabling Result %d, C3 Energy\n\n",j)
    end
end

end