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

% Z00UncAllLaunchEGA = (data.AsteroidUncertainty == 0 & string(data.OutTrajType) == "EGA");
% N01UncAllLaunchEGA = (data.AsteroidUncertainty == -0.1 & string(data.OutTrajType) == "EGA");
% P01UncAllLaunchEGA = (data.AsteroidUncertainty == 0.1 & string(data.OutTrajType) == "EGA");
% 
% Z00UncAllLaunchDIR = (data.AsteroidUncertainty == 0 & string(data.OutTrajType) == "Direct");
% N01UncAllLaunchDIR = (data.AsteroidUncertainty == -0.1 & string(data.OutTrajType) == "Direct");
% P01UncAllLaunchDIR = (data.AsteroidUncertainty == 0.1 & string(data.OutTrajType) == "Direct");
% figure(1)
% % plot(data.LaunchDate(Z00UncAllLaunch),data.OutFinDV(Z00UncAllLaunch),'bo')
% hold on
% % plot(data.LaunchDate(Z00UncAllLaunch),data.OutImpDV(Z00UncAllLaunch),'ro')
% plot(data.LaunchDate(Z00UncAllLaunchEGA),data.OutTotDV(Z00UncAllLaunchEGA),'bo')
% plot(data.LaunchDate(Z00UncAllLaunchDIR),data.OutTotDV(Z00UncAllLaunchDIR),'ro')
% xlabel('Launch Date')
% ylabel('\DeltaV (km/s)')
% title('Launch Date vs. \DeltaV Ast. Unc. = 0.00')
% legend('EGA Outbound','Direct Outbound','Location','best')
% 
% figure(2)
% % plot(data.LaunchDate(N01UncAllLaunchEGA),data.OutFinDV(N01UncAllLaunchEGA),'bo')
% hold on
% % plot(data.LaunchDate(N01UncAllLaunchEGA),data.OutImpDV(N01UncAllLaunchEGA),'ro')
% plot(data.LaunchDate(N01UncAllLaunchEGA),data.OutTotDV(N01UncAllLaunchEGA),'go')
% plot(data.LaunchDate(N01UncAllLaunchDIR),data.OutTotDV(N01UncAllLaunchDIR),'ro')
% xlabel('Launch Date')
% ylabel('\DeltaV (km/s)')
% title('Launch Date vs. \DeltaV Ast. Unc. = -0.1')
% legend('EGA Outbound','Direct Outbound','Location','best')
% 
% figure(3)
% % plot(data.LaunchDate(P01UncAllLaunchEGA),data.OutFinDV(P01UncAllLaunchEGA),'bo')
% hold on
% % plot(data.LaunchDate(P01UncAllLaunchEGA),data.OutImpDV(P01UncAllLaunchEGA),'ro')
% plot(data.LaunchDate(P01UncAllLaunchEGA),data.OutTotDV(P01UncAllLaunchEGA),'go')
% plot(data.LaunchDate(P01UncAllLaunchDIR),data.OutTotDV(P01UncAllLaunchDIR),'ro')
% xlabel('Launch Date')
% ylabel('\DeltaV (km/s)')
% title('Launch Date vs. \DeltaV Ast. Unc. = 0.1')
% legend('EGA Outbound','Direct Outbound','Location','best')

figure(1)
EGADates = unique(data.EgaDate);    
ImpDVsEGA= zeros(length(EGADates),3);    
FinDVsEGA = zeros(length(EGADates),3);
TotDVsEGA = zeros(length(EGADates),3);

for j = 1:length(EGADates)

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == 0     & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ))
        ImpDVsEGA(j,1) = max(data.OutImpDV(data.AsteroidUncertainty == 0   & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  )); 
    else
        ImpDVsEGA(j,1) = NaN;
    end

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == 0.1    & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ))
        ImpDVsEGA(j,2) = max(data.OutImpDV(data.AsteroidUncertainty == 0.1     & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ));  
    else
        ImpDVsEGA(j,2) = NaN;
    end

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == -0.1   & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ))
        ImpDVsEGA(j,3) = max(data.OutImpDV(data.AsteroidUncertainty == -0.1    & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ));  
    else
        ImpDVsEGA(j,3) = NaN;
    end
    
    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == 0     & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ))
        FinDVsEGA(j,1) = max(data.OutFinDV(data.AsteroidUncertainty == 0   & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  )); 
    else
        FinDVsEGA(j,1) = NaN;
    end

    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == 0.1    & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ))
        FinDVsEGA(j,2) = max(data.OutFinDV(data.AsteroidUncertainty == 0.1     & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ));  
    else
        FinDVsEGA(j,2) = NaN;
    end

    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == -0.1   & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ))
        FinDVsEGA(j,3) = max(data.OutFinDV(data.AsteroidUncertainty == -0.1    & string(data.OutTrajType) == 'EGA' & data.EgaDate == EGADates(j)  ));  
    else
        FinDVsEGA(j,3) = NaN;
    end
    
    ImpDVsEGA(j,:) = sort(ImpDVsEGA(j,:));
    FinDVsEGA(j,:) = sort(FinDVsEGA(j,:));
    TotDVsEGA(j,:) = ImpDVsEGA(j,:)+FinDVsEGA(j,:);
end


LaunchDates = unique(data.LaunchDate);    
ImpDVsDIR = zeros(length(LaunchDates),3);    
FinDVsDIR = zeros(length(LaunchDates),3);
TotDVsDIR = zeros(length(LaunchDates),3);

for j = 1:length(LaunchDates)

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == 0     & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        ImpDVsDIR(j,1) = max(data.OutImpDV(data.AsteroidUncertainty == 0   & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  )); 
    else
        ImpDVsDIR(j,1) = NaN;
    end

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == 0.1    & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        ImpDVsDIR(j,2) = max(data.OutImpDV(data.AsteroidUncertainty == 0.1     & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ));  
    else
        ImpDVsDIR(j,2) = NaN;
    end

    if ~isempty(data.OutImpDV(data.AsteroidUncertainty == -0.1   & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        ImpDVsDIR(j,3) = max(data.OutImpDV(data.AsteroidUncertainty == -0.1    & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ));  
    else
        ImpDVsDIR(j,3) = NaN;
    end
    
    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == 0     & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        FinDVsDIR(j,1) = max(data.OutFinDV(data.AsteroidUncertainty == 0   & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  )); 
    else
        FinDVsDIR(j,1) = NaN;
    end

    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == 0.1    & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        FinDVsDIR(j,2) = max(data.OutFinDV(data.AsteroidUncertainty == 0.1     & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ));  
    else
        FinDVsDIR(j,2) = NaN;
    end

    if ~isempty(data.OutFinDV(data.AsteroidUncertainty == -0.1   & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ))
        FinDVsDIR(j,3) = max(data.OutFinDV(data.AsteroidUncertainty == -0.1    & string(data.OutTrajType) == 'Direct' & data.LaunchDate == LaunchDates(j)  ));  
    else
        FinDVsDIR(j,3) = NaN;
    end
    
    ImpDVsDIR(j,:) = sort(ImpDVsDIR(j,:));
    FinDVsDIR(j,:) = sort(FinDVsDIR(j,:));
    TotDVsDIR(j,:) = ImpDVsDIR(j,:)+FinDVsDIR(j,:);
end



hold on
grid on

ynegImpDIR = abs(ImpDVsDIR(:,2)-ImpDVsDIR(:,1));
yposImpDIR = abs(ImpDVsDIR(:,3)-ImpDVsDIR(:,2));
ynegFinDIR = abs(FinDVsDIR(:,2)-FinDVsDIR(:,1));
yposFinDIR = abs(FinDVsDIR(:,3)-FinDVsDIR(:,2));
ynegTotDIR = ynegImpDIR + ynegFinDIR;
yposTotDIR = yposImpDIR + yposFinDIR;

ynegImpEGA = abs(ImpDVsEGA(:,2)-ImpDVsEGA(:,1));
yposImpEGA = abs(ImpDVsEGA(:,3)-ImpDVsEGA(:,2));
ynegFinEGA = abs(FinDVsEGA(:,2)-FinDVsEGA(:,1));
yposFinEGA = abs(FinDVsEGA(:,3)-FinDVsEGA(:,2));
ynegTotEGA = ynegImpEGA + ynegFinEGA;
yposTotEGA = yposImpEGA + yposFinEGA;

xneg = [];
xpos = [];

errorbar(datenum(LaunchDates),TotDVsDIR(:,2),ynegTotDIR,yposTotDIR,xneg,xpos,'ro')
errorbar(datenum(EGADates),TotDVsEGA(:,2),ynegTotEGA,yposTotEGA,xneg,xpos,'go')
datetick('x','mmm yyyy', 'keepticks', 'keeplimits')
xlabel('Launch Date')
ylabel('\DeltaV (km/s)')
title('Launch Date vs. \DeltaV ')
legend('Direct','EGA')




