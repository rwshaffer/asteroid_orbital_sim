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

Z00UncAllLaunch = (data.AsteroidUncertainty == 0 );
N01UncAllLaunch = (data.AsteroidUncertainty == -0.1 );
P01UncAllLaunch = (data.AsteroidUncertainty == 0.1 );
% 
% figure(1)
% plot(data.LaunchDate(Z00UncAllLaunch),data.OutFinDV(Z00UncAllLaunch),'bo')
% hold on
% plot(data.LaunchDate(Z00UncAllLaunch),data.OutImpDV(Z00UncAllLaunch),'ro')
% plot(data.LaunchDate(Z00UncAllLaunch),data.OutTotDV(Z00UncAllLaunch),'go')
% xlabel('Launch Date')
% ylabel('\DeltaV (km/s)')
% title('Launch Date vs. \DeltaV Ast. Unc. = 0.00')
% 
% figure(2)
% plot(data.LaunchDate(N01UncAllLaunch),data.OutFinDV(N01UncAllLaunch),'bo')
% hold on
% plot(data.LaunchDate(N01UncAllLaunch),data.OutImpDV(N01UncAllLaunch),'ro')
% plot(data.LaunchDate(N01UncAllLaunch),data.OutTotDV(N01UncAllLaunch),'go')
% xlabel('Launch Date')
% ylabel('\DeltaV (km/s)')
% title('Launch Date vs. \DeltaV Ast. Unc. = -0.1')
% 
% figure(3)
% plot(data.LaunchDate(P01UncAllLaunch),data.OutFinDV(P01UncAllLaunch),'bo')
% hold on
% plot(data.LaunchDate(P01UncAllLaunch),data.OutImpDV(P01UncAllLaunch),'ro')
% plot(data.LaunchDate(P01UncAllLaunch),data.OutTotDV(P01UncAllLaunch),'go')
% xlabel('Launch Date')
% ylabel('\DeltaV (km/s)')
% title('Launch Date vs. \DeltaV Ast. Unc. = 0.1')


figure(1)
LaunchDates = unique(data.LaunchDate);    
ImpDVs = zeros(length(LaunchDates),3);    
FinDVs = zeros(length(LaunchDates),3);
TotDVs = zeros(length(LaunchDates),3);

for j = 1:length(LaunchDates)

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == 0     & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        ImpDVs(j,1) = max(data.OutImpDV(data.AsteroidUncertainty == 0   & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  )); 
    else
        ImpDVs(j,1) = NaN;
    end

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == 0.1    & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        ImpDVs(j,2) = max(data.OutImpDV(data.AsteroidUncertainty == 0.1     & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ));  
    else
        ImpDVs(j,2) = NaN;
    end

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == -0.1   & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        ImpDVs(j,3) = max(data.OutImpDV(data.AsteroidUncertainty == -0.1    & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ));  
    else
        ImpDVs(j,3) = NaN;
    end
    
    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == 0     & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        FinDVs(j,1) = max(data.OutFinDV(data.AsteroidUncertainty == 0   & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  )); 
    else
        FinDVs(j,1) = NaN;
    end

    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == 0.1    & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        FinDVs(j,2) = max(data.OutFinDV(data.AsteroidUncertainty == 0.1     & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ));  
    else
        FinDVs(j,2) = NaN;
    end

    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == -0.1   & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        FinDVs(j,3) = max(data.OutFinDV(data.AsteroidUncertainty == -0.1    & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ));  
    else
        FinDVs(j,3) = NaN;
    end
    
    ImpDVs(j,:) = sort(ImpDVs(j,:));
    FinDVs(j,:) = sort(FinDVs(j,:));
    TotDVs(j,:) = ImpDVs(j,:)+FinDVs(j,:);
end



hold on
grid on
ynegImp = abs(ImpDVs(:,2)-ImpDVs(:,1));
yposImp = abs(ImpDVs(:,3)-ImpDVs(:,2));
ynegFin = abs(FinDVs(:,2)-FinDVs(:,1));
yposFin = abs(FinDVs(:,3)-FinDVs(:,2));
ynegTot = ynegImp + ynegFin;
yposTot = yposImp + yposFin;
xneg = [];
xpos = [];
errorbar(datenum(LaunchDates),ImpDVs(:,2),ynegImp,yposImp,xneg,xpos,'o')
errorbar(datenum(LaunchDates),FinDVs(:,2),ynegFin,yposFin,xneg,xpos,'o')
errorbar(datenum(LaunchDates),TotDVs(:,2),ynegTot,yposTot,xneg,xpos,'o')
datetick('x','dd-mmm-YYYY')
xlabel('Launch Date')
ylabel('\DeltaV (km/s)')
title('Launch Date vs. \DeltaV - Direct Outbound')
legend('Electric \DeltaV','Impulsive \DeltaV','Total Outbound \DeltaV','Location','best')

