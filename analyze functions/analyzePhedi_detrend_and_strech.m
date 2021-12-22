function analyzePhedi_detrend_and_strech(DataStruct, varargin)
% PhediAnalyze_detrend_and_strech will plot normalized data in different
% configurations as requested by user.
%
% varargin can contain:
%   a vector signifying the location to consider in this presentation.
%   name of figures to present. These will include (CASE SENSITIVE):
%       'TrendFits'
%       'TrendNormSmth'
%       'RawNormSmth'
%       'RawNorm'
%       'RawNormMesh'


VararginVectors = cellfun(@isnumeric, varargin);
plotNames = varargin(~VararginVectors);

PhediLocation = DataStruct.PhediData.PhediLocation;
timeVec = DataStruct.PhediData.timeVec;

if nnz(VararginVectors)>0
    boundaries = varargin{find(VararginVectors,1,'first')};
    PhediLocation = PhediLocation(boundaries,:);
    timeVec = timeVec(boundaries);
end
PhediLocation_omitnan = PhediLocation(~isnan(PhediLocation));
FigColors = MyVaryColor(size(PhediLocation,2));

P = [];
PhediLocation_detrend = [];
PhediLocation_streched = [];
PhediLocation_smoothed = [];
PhediLocation_raw_stretch = [];
PhediLocation_raw_smth = [];
for i=1:size(PhediLocation,2)
    %--- create no nan logical
    noNanLogic = ~isnan(PhediLocation(:,i));
    %--- fit a polynomial
    p = polyfit(timeVec(noNanLogic),PhediLocation(noNanLogic,i),2);
    P = cat(2,P,p(:));
    %--- remove linear trend
    LinFit = polyval(p,timeVec);
    PhediLocation_detrend_i = PhediLocation(:,i)-LinFit(:);
    %--- strech between -1 to 1
    PhediLocation_dtrnd_streched_i = stretchAmpMinus1to1(PhediLocation_detrend_i);
    %--- smooth signal
    smoothes_streched_i = supsmu(timeVec(noNanLogic),PhediLocation_dtrnd_streched_i(noNanLogic));
    %--- stretch raw
    PhediLocation_raw_streched_i = stretchAmpMinus1to1(PhediLocation(:,i));
    %--- smoothed streched raw
    PhediLocation_raw_smth_i = supsmu(timeVec(noNanLogic),PhediLocation_raw_streched_i(noNanLogic));
    
    %--- save the data to matrixes
    PhediLocation_detrend = cat(2,PhediLocation_detrend,PhediLocation_detrend_i);
    PhediLocation_streched = cat(2,PhediLocation_streched,PhediLocation_dtrnd_streched_i);
    PhediLocation_smoothed = cat(2,PhediLocation_smoothed,smoothes_streched_i);
    PhediLocation_raw_stretch = cat(2,PhediLocation_raw_stretch,PhediLocation_raw_streched_i);
    PhediLocation_raw_smth = cat(2,PhediLocation_raw_smth,PhediLocation_raw_smth_i);
    
    
end

%% plot stuff


% 'TrendFits'
% 'TrendNormSmth'
% 'RawNormSmth'
% 'RawNorm'
% 'RawNormMesh'

%--- create figure for polynom presenting
if nnz(strcmpi(plotNames,'TrendFits'))>0
    polyFig = figure;
    hold on;
    for i=1:size(PhediLocation,2)
        if DataStruct.PhediData.slopeIncline(i)>0
            plot(timeVec,PhediLocation(:,i),'*-','Color',FigColors(i,:));
            plot(timeVec,polyval(P(:,i),timeVec),'LineWidth',1.5,'Color',FigColors(i,:));
        elseif DataStruct.PhediData.slopeIncline(i)<0
            plot(timeVec,PhediLocation(:,i),'.-','Color',FigColors(i,:));
            plot(timeVec,polyval(P(:,i),timeVec),'LineWidth',1.5,'Color',FigColors(i,:));
        else
            plot(timeVec,PhediLocation(:,i),'--','Color',FigColors(i,:));
            plot(timeVec,polyval(P(:,i),timeVec),'LineWidth',1.5,'Color',FigColors(i,:));
        end
    end
    title([DataStruct.BigPicRotStruct.details, 'fits for trending']);
    xlabel('t [s]');
    ylabel('location [m]');
end

%--- create figure for detrend
if nnz(strcmpi(plotNames,'TrendNormSmth'))>0
    dtrndFig = figure;
    hold on;
    for i=1:size(PhediLocation,2)
        if DataStruct.PhediData.slopeIncline(i)>0
            plot(timeVec(noNanLogic),PhediLocation_smoothed(:,i),'*-','Color',FigColors(i,:));
        elseif DataStruct.PhediData.slopeIncline(i)<0
            plot(timeVec(noNanLogic),PhediLocation_smoothed(:,i),'.-','Color',FigColors(i,:));
        else
            plot(timeVec(noNanLogic),PhediLocation_smoothed(:,i),'LineWidth',1.5,'Color',FigColors(i,:));
        end
    end
    title({'trended and smoothed',DataStruct.BigPicRotStruct.details});
    xlabel('t [s]');
    ylabel('\Deltalocation [m]');
end

%--- create figure for smoothed raw
if nnz(strcmpi(plotNames,'RawNormSmth'))>0
    rawSmthdFig = figure;
    hold on;
    for i=1:size(PhediLocation,2)
        if DataStruct.PhediData.slopeIncline(i)>0
            plot(timeVec(noNanLogic),PhediLocation_raw_smth(:,i),'*-','Color',FigColors(i,:));
        elseif DataStruct.PhediData.slopeIncline(i)<0
            plot(timeVec(noNanLogic),PhediLocation_raw_smth(:,i),'.-','Color',FigColors(i,:));
        else
            plot(timeVec(noNanLogic),PhediLocation_raw_smth(:,i),'LineWidth',1.5,'Color',FigColors(i,:));
        end
    end
    title({'raw data smoothed and normalized', DataStruct.BigPicRotStruct.details});
    xlabel('t [s]');
    ylabel('amp [arb]');
end

%--- create figure for raw streched
if nnz(strcmpi(plotNames,'RawNorm'))>0
    rawStrchFig = figure;
    hold on;
    for i=1:size(PhediLocation,2)
        if DataStruct.PhediData.slopeIncline(i)>0
            plot(timeVec,PhediLocation_raw_stretch(:,i),'*-','Color',FigColors(i,:));
        elseif DataStruct.PhediData.slopeIncline(i)<0
            plot(timeVec,PhediLocation_raw_stretch(:,i),'.-','Color',FigColors(i,:));
        else
            plot(timeVec,PhediLocation_raw_stretch(:,i),'LineWidth',1.5,'Color',FigColors(i,:));
        end
    end
    title({'raw data normalized', DataStruct.BigPicRotStruct.details});
    xlabel('t [s]');
    ylabel('amp [arb]');
end

%--- plot mesh of normalized raw data
if nnz(strcmpi(plotNames,'RawNormMesh'))>0
    slopeLogical = DataStruct.PhediData.slopeIncline>0;
    [X,Y] = meshgrid(timeVec,1:size(PhediLocation_raw_stretch(:,slopeLogical),2));
    FN = DataStruct.SgData.N;
    FS = DataStruct.SgData.F;
    ExpState = ['Nmean=',num2str(round(mean(FN))),'N1=',num2str(round(mean(FN(1:100)))),' Nend=',num2str(round(mean(FN(end-100:end)))),...
        '  Fmean=',num2str(round(mean(FS))),'F1=',num2str(round(mean(FS(1:100)))),' Fend=',num2str(round(mean(FS(end-100:end))))];
    figure;
    mesh(X,Y,PhediLocation_raw_stretch(:,slopeLogical)');
    colormap(jet);
    title({'mesh of raw data normalized', DataStruct.BigPicRotStruct.details,ExpState});
    xlabel('t [s]');
    ylabel('position on block [arb]');
    zlabel('amp [arb]');
    
end

end