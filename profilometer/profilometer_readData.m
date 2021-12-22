function [fixed_xx, fixedAlt, fixedRowNums, ScanParameters] = profilometer_readData(folderPath, sampNum, ifPlot)
%profilometer_readData reads data of profilometer
% sampNum is the number of measurment in the folder specified, when
% organized by accending name.


if nargin<3
    ifPlot = 0;
end

%% read data
test=dir_names(folderPath);
test = test(1:2:end);
MeasureName=test{sampNum}(1:end-4);

DataFile = [MeasureName,'.txt'];
ParameterFile = [MeasureName,'_parameters.txt'];

ImportedData=importdata([folderPath,'\', DataFile]);

RowNums = ImportedData.data(:,5);
RawAlt=ImportedData.data(:,2);
xx=ImportedData.data(:,3);

%% read parameters
ScanParameters = profilometer_readParameters([folderPath,'\', ParameterFile]);

%% fix data
noiseThresh = 59;
fixedAlt = RawAlt(RawAlt>noiseThresh);
fixed_xx = xx(RawAlt>noiseThresh);
fixedRowNums = RowNums(RawAlt>noiseThresh);

%--- number of rows scanned:
rowsQuant = length(unique(fixedRowNums));
if rowsQuant~=1
    multipleRowFlag = 1;
else
    multipleRowFlag = 0;
end

%--- linear remove trend:
% fixed_xx = xx;
% fixedAlt = RawAlt;
% fixed_Altitued = detrend(fixed_Altitued);

%--- remove trend polynomial trend:
polDeg = 1;
smoothPar = 1000;
for rowIdx = 1:rowsQuant
    %--- work on current row:
    currentAlt = fixedAlt(fixedRowNums==rowIdx);
    currentXX = fixed_xx(fixedRowNums==rowIdx);
    
    %--- fix:
    smoothedAlt= smooth(currentAlt,smoothPar);
    p = polyfit(currentXX,smoothedAlt,polDeg);
    currentAlt = currentAlt-polyval(p,currentXX);
    
    %--- chnge the data:
    fixedAlt(fixedRowNums==rowIdx)=-currentAlt;
end

%% plot
if ifPlot
    figure; hold all
    
    if multipleRowFlag
        %         uniqLocs = unique(fixedRowNums);
        LGND = {};
        for rowIdx = 1:rowsQuant
            plot(fixed_xx(fixedRowNums==rowIdx)', fixedAlt(fixedRowNums==rowIdx)');
            LGND = cat(1,LGND,num2str(rowIdx*ScanParameters.RowStep.Val));
        end
        legend(LGND);
    else
        plot(xx,RawAlt,'.')
        plot(fixed_xx-fixed_xx(1),fixedAlt)
        plot(fixed_xx, polyval(p,fixed_xx));
        legend({'raw data', 'fixed data'});
    end
    title(MeasureName(1:end-20));
    xlabel('x [\mum]')
    ylabel('Altitude [\mum]')
end

end