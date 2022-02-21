function phase = time_to_phase(time

% INPUTS
% time - date as string in most date formats, use mm/dd/yyyy as default
%
% OUTPUTS
% phase - phase angle

phasedata = readtable('EarthTargetPhaseAngle.csv');
phasedata.Time = datenum(phasedata.Time);
timeNum = datenum(time);

for i = 1:length(phasedata.Time)-1
    if phasedata.Time(i)<= timeNum & timeNum < phasedata.Time(i+1)
        phase = interp1([phasedata.Time(i),phasedata.Time(i+1)],[phasedata.Phase(i),phasedata.Phase(i+1)],timeNum);
        break
    end
end
end