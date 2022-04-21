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
data.LaunchDate = datetime(data.LaunchDate,'InputFormat','dd MMM yyyy HH:mm:ss','Format','dd/MM/yyyy');
data.EgaDate = datetime(data.EgaDate,'InputFormat','dd MMM yyyy HH:mm:ss','Format','dd/MM/yyyy');
data.AsteroidDepartureDate = datetime(data.AsteroidDepartureDate,'InputFormat','dd MMM yyyy HH:mm:ss','Format','dd/MM/yyyy');
data.EarthArrivalDate = datetime(data.EarthArrivalDate,'InputFormat','dd MMM yyyy HH:mm:ss','Format','dd/MM/yyyy');

data.AsteroidDepartureDate = dateshift(data.AsteroidDepartureDate,'start','week');

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

Z00UncAllLaunch = (data.AsteroidUncertainty == 0 );
N01UncAllLaunch = (data.AsteroidUncertainty == -0.1 );
P01UncAllLaunch = (data.AsteroidUncertainty == 0.1 );
% 
% figure(1)
% plot(data.AsteroidDepartureDate(Z00UncAllLaunch),data.RetDV(Z00UncAllLaunch),'bo')
% hold on
% plot(data.AsteroidDepartureDate(Z00UncAllLaunch),data.TotDV(Z00UncAllLaunch),'ro')
% xlabel('Asteroid Departure Date')
% ylabel('\DeltaV (km/s)')
% title('Asteroid Departure Date vs. \DeltaV Ast. Unc. = 0.00')
% 
% figure(2)
% plot(data.AsteroidDepartureDate(N01UncAllLaunch),data.RetDV(N01UncAllLaunch),'bo')
% hold on
% plot(data.AsteroidDepartureDate(N01UncAllLaunch),data.TotDV(N01UncAllLaunch),'ro')
% xlabel('Asteroid Departure Date')
% ylabel('\DeltaV (km/s)')
% title('Asteroid Departure Date vs. \DeltaV Ast. Unc. = -0.1')
% 
% figure(3)
% plot(data.AsteroidDepartureDate(P01UncAllLaunch),data.RetDV(P01UncAllLaunch),'bo')
% hold on
% plot(data.AsteroidDepartureDate(P01UncAllLaunch),data.TotDV(P01UncAllLaunch),'ro')
% xlabel('Asteroid Departure Date')
% ylabel('\DeltaV (km/s)')
% title('Asteroid Departure Date vs. \DeltaV Ast. Unc. = 0.1')


figure(1)
DepDates = data.AsteroidDepartureDate;   
RetDVs = zeros(length(DepDates),3);    


for j = 1:length(DepDates)

    if ~isempty(data.RetDV(data.AsteroidUncertainty == 0     & data.AsteroidDepartureDate == DepDates(j)  ))
        RetDVs(j,1) = max(data.RetDV(data.AsteroidUncertainty == 0   & data.AsteroidDepartureDate == DepDates(j)  )); 
    else
        RetDVs(j,1) = NaN;
    end

    if ~isempty(data.RetDV(data.AsteroidUncertainty == 0.1    & data.AsteroidDepartureDate == DepDates(j)  ))
        RetDVs(j,2) = max(data.RetDV(data.AsteroidUncertainty == 0.1     & data.AsteroidDepartureDate == DepDates(j)  ));  
    else
        RetDVs(j,2) = NaN;
    end

    if ~isempty(data.RetDV(data.AsteroidUncertainty == -0.1   & data.AsteroidDepartureDate == DepDates(j)  ))
        RetDVs(j,3) = max(data.RetDV(data.AsteroidUncertainty == -0.1    & data.AsteroidDepartureDate == DepDates(j)  ));  
    else
        RetDVs(j,3) = NaN;
    end
    
        
    RetDVs(j,:) = sort(RetDVs(j,:));
end



hold on
grid on
yneg = abs(RetDVs(:,2)-RetDVs(:,1));
ypos = abs(RetDVs(:,3)-RetDVs(:,2));
xneg = [];
xpos = [];
errorbar(datenum(DepDates),RetDVs(:,1),yneg,ypos,xneg,xpos,'o')
xlim([7.403e5 7.411e5])
datetick('x','dd-mmm-YYYY','keeplimits')
xlabel('Departure Date')
ylabel('\DeltaV (km/s)')
title('Departure Date vs. \DeltaV - Direct Return')
