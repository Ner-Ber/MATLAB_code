% risngFig = figure; hold on;
% yldTauFig = figure; hold on;
% BigSlpFig = figure; hold on;
% yldTauFigByKink = figure; hold on;

D = setdiff(dir_names,my_dir);
Smega = struct([]);
for d = 1:length(D)
    load(D{d});
    relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
    relevantEvents = relevantEvents(relevantEvents~=1);
    try
%         risingX = nan(length(relevantEvents),1);
%         yldTvec1st = nan(length(relevantEvents),1);
%         BigSlpVec = nan(length(relevantEvents),1);
        yldTvec2nd = nan(length(relevantEvents),1);
        Cf_vec = nan(length(relevantEvents),1);
        for I=1:length(relevantEvents)
            try
                DataStruct = PhediStructCell{relevantEvents(I)};
                DataStruct.PhediData = PhediStructCell{relevantEvents(I)}.PhediDataSimpSmth;
                Cf_vec(I) = DataStruct.PhediData.Cf;
                AvgPhediStruct  = phedi_averagePhedi(DataStruct);
                S = yield_findYieldStressStrain(DataStruct);
                S.Date = PhediStructCell{relevantEvents(1)}.ExperimentData.ExpDate;
                S.Hour = PhediStructCell{relevantEvents(1)}.ExperimentData.ExpHour;
                S.Event = relevantEvents(I);
                Smega = cat(1,Smega,S);
%                 risingX(I) = S.risingVelPoint;
%                 yldTvec(I) = S.yieldSxy1st;
%                 BigSlpVec(I) = S.scndSlope;
%                 yldTvec2nd(I) = S.yieldSxy2nd;
            catch
            end
        end
%         figure(risngFig);
%         plot(Cf_vec,risingX,'*','DisplayName',[PhediStructCell{relevantEvents(1)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(1)}.ExperimentData.ExpHour]);
%         figure(yldTauFig);
%         plot(Cf_vec,yldTvec1st,'*','DisplayName',[PhediStructCell{relevantEvents(1)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(1)}.ExperimentData.ExpHour]);
%         figure(BigSlpFig);
%         plot(Cf_vec,BigSlpVec,'*','DisplayName',[PhediStructCell{relevantEvents(1)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(1)}.ExperimentData.ExpHour]);
%         figure(yldTauFigByKink);
%         plot(Cf_vec,yldTvec2nd,'*','DisplayName',[PhediStructCell{relevantEvents(1)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(1)}.ExperimentData.ExpHour]);
%         pause(3);
    catch
    end
end
