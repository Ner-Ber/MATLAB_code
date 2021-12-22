StructCellList = ...
    {'PhediStructCell_1to26_Exp1','PhediStructCell_1to24_Exp2','PhediStructCell_1to30_Exp3','PhediStructures_1to30_Exp4','PhediStructures_1to28_Exp5'};
markerList = {'o','+','*','x','s','d','^','v','>','<','p','h'};
colorList = {rgb('Blue'),rgb('Red'),rgb('Orange'),rgb('Purple'),rgb('Green'),...
    rgb('HotPink'), rgb('Yellow'),rgb('Gray')};

h = zeros(size(StructCellList));
NameN = {};
figure;
hold on;
for s = 1:length(StructCellList)
    %% prepare data
    PhediStructCell = eval(StructCellList{s});
    AmpCell = cell(length(PhediStructCell),1);
    F_vec = zeros(length(PhediStructCell),1);
    N_vec = zeros(length(PhediStructCell),1);
    AmpMat = nan(40,length(PhediStructCell));
    
    for i = 1:length(PhediStructCell)
        try
            PhediLocations = PhediStructCell{i}.PhediData.PhediLocation;
            timeVec = PhediStructCell{i}.PhediData.timeVec;
            logicalTime = timeVec<-0.2e-3;
            
            PhediLocations = PhediLocations(logicalTime,:);
            timeVec = timeVec(logicalTime);
            
            PhediLocationAmp = max(PhediLocations,[],1,'omitnan') - min(PhediLocations,[],1,'omitnan');
            AmpCell{i} = PhediLocationAmp;
            AmpMat(1:length(PhediLocationAmp),i) = PhediLocationAmp(:);
            F_vec(i) = mean(PhediStructCell{i}.SgData.F);
            N_vec(i) = mean(PhediStructCell{i}.SgData.N);
            
            %         plot(F_vec(i),PhediLocationAmp,'Color',colorList{s},'Marker',markerList{i});
            %         boxplot(PhediLocationAmp,F_vec(i),'Colors',colorList{s},'PlotStyle','compact','Whisker',inf);
        catch
        end
    end
    %     Pos = round(repmat(F_vec(:)',size(AmpMat,1),1));
    %     Pos = [Pos(:); 300;-5];
    %     boxplot([AmpMat(:); nan; nan],Pos,'Colors',colorList{s},'PlotStyle','compact','Whisker',inf,'Positions',Pos);
    
    meanVals = mean(abs(AmpMat),'omitnan');
    upperBound = abs(meanVals-max(abs(AmpMat),[],1,'omitnan'));
    lowerBound = abs(meanVals-min(abs(AmpMat),[],1,'omitnan'));
    errorbar(F_vec,meanVals,lowerBound,upperBound,'*','Color',colorList{s});
    NameN{s} = ['N=',num2str(mean(N_vec(3:end)))];
end

legend(NameN);