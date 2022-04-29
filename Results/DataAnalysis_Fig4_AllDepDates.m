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
data.LaunchDate = datetime(data.LaunchDate,'InputFormat','dd MMM yyyy HH:mm:ss');
data.EgaDate = datetime(data.EgaDate,'InputFormat','dd MMM yyyy HH:mm:ss');
data.AsteroidDepartureDate = datetime(data.AsteroidDepartureDate,'InputFormat','dd MMM yyyy HH:mm:ss');
data.EarthArrivalDate = datetime(data.EarthArrivalDate,'InputFormat','dd MMM yyyy HH:mm:ss');

data.AsteroidDepartureDate = dateshift(data.AsteroidDepartureDate,'start','day');

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

[data.AsteroidDepartureDate,sortIdx] = sort(data.AsteroidDepartureDate);
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
data.LaunchDate = data.LaunchDate(sortIdx);
data.EarthArrivalDate = data.EarthArrivalDate(sortIdx);

data.TotTOF = data.RetTOF+data.OutTOF;

Z00UncAllLaunch = (data.AsteroidUncertainty == 0 );
N01UncAllLaunch = (data.AsteroidUncertainty == -0.1 );
P01UncAllLaunch = (data.AsteroidUncertainty == 0.1 );

% figure(1)
% % plot(data.AsteroidDepartureDate(Z00UncAllLaunch),data.OutTOF(Z00UncAllLaunch),'bo')
% hold on
% plot(data.AsteroidDepartureDate(Z00UncAllLaunch),data.RetTOF(Z00UncAllLaunch),'ro')
% % plot(data.AsteroidDepartureDate(Z00UncAllLaunch),data.TotTOF(Z00UncAllLaunch),'go')
% xlabel('Asteroid Departure Date')
% ylabel('ToF (days)')
% title('Asteroid Departure Date vs. TOF Ast. Unc. = 0.00')
% 
% figure(2)
% % plot(data.AsteroidDepartureDate(N01UncAllLaunch),data.OutTOF(N01UncAllLaunch),'bo')
% hold on
% plot(data.AsteroidDepartureDate(N01UncAllLaunch),data.RetTOF(N01UncAllLaunch),'ro')
% % plot(data.AsteroidDepartureDate(N01UncAllLaunch),data.TotTOF(N01UncAllLaunch),'go')
% xlabel('Asteroid Departure Date')
% ylabel('ToF (days)')
% title('Asteroid Departure Date vs. TOF Ast. Unc. = -0.1')
% 
% figure(3)
% % plot(data.AsteroidDepartureDate(P01UncAllLaunch),data.OutTOF(P01UncAllLaunch),'bo')
% hold on
% plot(data.AsteroidDepartureDate(P01UncAllLaunch),data.RetTOF(P01UncAllLaunch),'ro')
% % plot(data.AsteroidDepartureDate(P01UncAllLaunch),data.TotTOF(P01UncAllLaunch),'go')
% xlabel('Asteroid Departure Date')
% ylabel('ToF (days)')
% title('Asteroid Departure Date vs. TOF Ast. Unc. = 0.1')

%% new

DepDates = unique(data.AsteroidDepartureDate);    
RetTOFs = zeros(length(DepDates),3);    

for j = 1:length(DepDates)

    if ~isempty(data.RetTOF(data.AsteroidUncertainty == 0     &  data.AsteroidDepartureDate == DepDates(j)  ))
        RetTOFs(j,1) = max(data.RetTOF(data.AsteroidUncertainty == 0     &  data.AsteroidDepartureDate == DepDates(j)  )); 
    else
        RetTOFs(j,1) = NaN;
    end

    if ~isempty(data.RetTOF(data.AsteroidUncertainty == 0.1     &  data.AsteroidDepartureDate == DepDates(j)  ))
        RetTOFs(j,2) = max(data.RetTOF(data.AsteroidUncertainty == 0.1     &  data.AsteroidDepartureDate == DepDates(j)  ));  
    else
        RetTOFs(j,2) = NaN;
    end

    if ~isempty(data.RetTOF(data.AsteroidUncertainty == -0.1     &  data.AsteroidDepartureDate == DepDates(j)  ))
        RetTOFs(j,3) = max(data.RetTOF(data.AsteroidUncertainty == -0.1     &  data.AsteroidDepartureDate == DepDates(j)  ));  
    else
        RetTOFs(j,3) = NaN;
    end

    RetTOFs(j,:) = sort(RetTOFs(j,:));
end

hold on
grid on
yneg = abs(RetTOFs(:,2)-RetTOFs(:,1));
ypos = abs(RetTOFs(:,3)-RetTOFs(:,2));
xneg = [];
xpos = [];
errorbar(datenum(DepDates),RetTOFs(:,2),yneg,ypos,xneg,xpos,'o')

datetick('x','dd-mmm-YYYY')
xlabel('Asteroid Departure Date')
ylabel('ToF (days)')
title('Launch Date vs. Time of Flight')



