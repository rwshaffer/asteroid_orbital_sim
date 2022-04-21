clear all
close all
clc
% varNames = {'LaunchDate','MiningDuration','OutTrajType','AsterioidUncertainty','OutFinDV','OutImpDV','RetDV','MaxThrust','C3Energy','OutTOF','RetTOF','TotDV','AsteroidDepartureDate','EarthArrivalDate'};
% dataStartLine = 2;
% opts = delimitedTextImportOptions('VariableNames',varNames,'DataLines', dataStartLine)


%% Importing Data from CSV File 

filename = "ExpRes.csv";

varNames = {'LaunchDate','MiningDuration','OutTrajType',...
            'AsteroidUncertainty','OutFinDV','OutImpDV',...
            'RetDV','MaxThrust','C3Energy','EgaDate','OutTOF',...
            'RetTOF','TotDV','AsteroidDepartureDate','EarthArrivalDate'};

varTypes = {'char','double','char','double','double','double','double',...
            'double','double','char','double','double','double','char','char'};
delimiter = ',';
dataStartLine = 2;
extraColRule = 'ignore';

opts = delimitedTextImportOptions('VariableNames',varNames,...
                                'VariableTypes',varTypes,...
                                'Delimiter',delimiter,...
                                'DataLines', dataStartLine,...
                                'ExtraColumnsRule',extraColRule);

data = readtable(filename,opts);

% datatime variables
data.LaunchDate = datetime(data.LaunchDate,'InputFormat','dd MMM yyyy HH:mm:ss','Format','dd MMM yyyy');
data.EgaDate = datetime(data.EgaDate,'InputFormat','dd MMM yyyy HH:mm:ss','Format','dd MMM yyyy');
data.AsteroidDepartureDate = datetime(data.AsteroidDepartureDate,'InputFormat','dd MMM yyyy HH:mm:ss','Format','dd MMM yyyy');
data.EarthArrivalDate = datetime(data.EarthArrivalDate,'InputFormat','dd MMM yyyy HH:mm:ss','Format','dd MMM yyyy');

data.OutFinDV = data.OutFinDV + sqrt(data.C3Energy);
data.OutTotDV = data.OutFinDV + data.OutImpDV;

[numData,numVars]=size(data);

% [LaunchDates,sortIdx] = sort(data.LaunchDate);
% MiningDurations = data.MiningDuration(sortIdx);
% OutTrajTypes = data.OutTrajType(sortIdx);
% AsteroidUncertainties = data.AsteroidUncertainty(sortIdx);
% OutFinDVs = data.OutFinDV(sortIdx);
% OutImpDVs = data.OutImpDV(sortIdx);
% OutTotDVs = data.OutTotDV(sortIdx);
% RetDVs = data.RetDV(sortIdx);
% C3Energies = data.C3Energy(sortIdx);
% EgaDates = data.EgaDate(sortIdx);
% OutTOFs = data.OutTOF(sortIdx);
% RetTOFs = data.RetTOF(sortIdx);
% TotDVs = data.TotDV(sortIdx);
% AsteroidDepartureDates = data.AsteroidDepartureDate(sortIdx);
% EarthArrivalDates = data.EarthArrivalDate(sortIdx);

[data.LaunchDate,sortIdx] = sort(data.LaunchDate);
data.MiningDuration = data.MiningDuration(sortIdx);
data.OutTrajType = data.OutTrajType(sortIdx);
data.AsteroidUncertainty = data.AsteroidUncertainty(sortIdx);
data.OutFinDV = data.OutFinDV(sortIdx);
data.OutImpDV = data.OutImpDV(sortIdx);
data.OutTotDV = data.OutTotDV(sortIdx);
data.RetDV = data.RetDV(sortIdx);
data.C3Energy = data.C3Energy(sortIdx);
data.EgaDate = data.EgaDate(sortIdx);
data.OutTOF = data.OutTOF(sortIdx);
data.RetTOF = data.RetTOF(sortIdx);
data.TotDV = data.TotDV(sortIdx);
data.AsteroidDepartureDate = data.AsteroidDepartureDate(sortIdx);
data.EarthArrivalDate = data.EarthArrivalDate(sortIdx);

data.TotTOF = data.RetTOF+data.OutTOF;

mineDurs = unique(data.MiningDuration);

figure(1)
k = 1;
legendEntries = string([]);

for i = 1:length(mineDurs)
    mineLaunchDates = data.LaunchDate(data.MiningDuration == mineDurs(i));    
    mineTOFs = zeros(length(mineLaunchDates),3);    

    for j = 1:length(mineLaunchDates)
        
        if ~isempty(data.TotTOF(data.AsteroidUncertainty == 0    & data.MiningDuration == mineDurs(i) & data.LaunchDate >= datetime("01-Jul-2025 00:00:00") & data.LaunchDate <= datetime("31-Jul-2025 23:59:59") & data.LaunchDate == mineLaunchDates(j)  ))
            mineTOFs(j,1) = max(data.TotTOF(data.AsteroidUncertainty == 0    & data.MiningDuration == mineDurs(i) & data.LaunchDate >= datetime("01-Jul-2025 00:00:00") & data.LaunchDate <= datetime("31-Jul-2025 23:59:59") & data.LaunchDate == mineLaunchDates(j)  )); 
        else
            mineTOFs(j,1) = NaN;
        end
        
        if ~isempty(data.TotTOF(data.AsteroidUncertainty == 0.1    & data.MiningDuration == mineDurs(i) & data.LaunchDate >= datetime("01-Jul-2025 00:00:00") & data.LaunchDate <= datetime("31-Jul-2025 23:59:59") & data.LaunchDate == mineLaunchDates(j)  ))
            mineTOFs(j,2) = max(data.TotTOF(data.AsteroidUncertainty == 0.1    & data.MiningDuration == mineDurs(i) & data.LaunchDate >= datetime("01-Jul-2025 00:00:00") & data.LaunchDate <= datetime("31-Jul-2025 23:59:59") & data.LaunchDate == mineLaunchDates(j)  ));  
        else
            mineTOFs(j,2) = NaN;
        end
        
        if ~isempty(data.TotTOF(data.AsteroidUncertainty == -0.1    & data.MiningDuration == mineDurs(i) & data.LaunchDate >= datetime("01-Jul-2025 00:00:00") & data.LaunchDate <= datetime("31-Jul-2025 23:59:59") & data.LaunchDate == mineLaunchDates(j)  ))
            mineTOFs(j,3) = max(data.TotTOF(data.AsteroidUncertainty == -0.1    & data.MiningDuration == mineDurs(i) & data.LaunchDate >= datetime("01-Jul-2025 00:00:00") & data.LaunchDate <= datetime("31-Jul-2025 23:59:59") & data.LaunchDate == mineLaunchDates(j)  ));  
        else
            mineTOFs(j,3) = NaN;
        end
        
        mineTOFs(j,:) = sort(mineTOFs(j,:));
    end
    
    legendEntries(k) = string(mineDurs(i))+" Days";
    k = k + 1;
    
    hold on
    grid on
    yneg = abs(mineTOFs(:,2)-mineTOFs(:,1));
    ypos = abs(mineTOFs(:,3)-mineTOFs(:,2));
    xneg = [];
    xpos = [];
    errorbar(datenum(mineLaunchDates),mineTOFs(:,2),yneg,ypos,xneg,xpos,'o')
end
datetick('x','dd-mmm-YYYY')
xlabel('Launch Date')
ylabel('ToF (days)')
title('Launch Date vs. Time of Flight - July Launch')
legend(legendEntries,'Location','bestoutside');

