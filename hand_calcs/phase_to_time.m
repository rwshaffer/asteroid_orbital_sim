function time = phase_to_time(phase)

% INPUTS
% phase - phase angle of enquiry
%
% OUTPUTS
% time - vector of dates at phase angle as strings

phasedata = readtable('EarthTargetPhaseAngle.csv');
phasedata.Time = datenum(phasedata.Time);

time = [];

for i = 1:length(phasedata.Time)-1
    if phasedata.Phase(i)<= phase & phase < phasedata.Phase(i+1)
        timeNum = interp1([phasedata.Phase(i),phasedata.Phase(i+1)],[phasedata.Time(i),phasedata.Time(i+1)],phase);
        time = [time;datestr(timeNum,'mm/dd/yyyy')];
    end
end
end