%%
% PhediStructCelld = [];
%**** load PhediStructCell****
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
FWHM_vec = [];
Cf_vec = [];
FWHM_Y_vec = [];
Loc2Right_vec = [];
Loc2Left_vec = [];
lengthbyLEFM = [];
for I=relevantEvents
    %-- change to different smoothing:
    DataStruct = phedi_smoothByfreqFilter(PhediStructCell{I});
    
    [FWHM,FWHM_Y,Loc2Right,Loc2Left] = phedi_FWHMofMean(DataStruct);
    FWHM_vec = cat(1,FWHM_vec,FWHM);
    FWHM_Y_vec = cat(1,FWHM_Y_vec,FWHM_Y);
    Cf_vec = cat(1,Cf_vec,DataStruct.PhediData.Cf);
    Loc2Right_vec = cat(1,Loc2Right_vec,Loc2Right);
    Loc2Left_vec = cat(1,Loc2Left_vec,Loc2Left);
    
    [lengthbyLEFM_i,~] = intersections(DataStruct.solAtInter.x,DataStruct.solAtInter.vx,[-100 100],2*FWHM_Y*[1 1]);
    if isempty(lengthbyLEFM_i)
        lengthbyLEFM_i = 0;
    end
    lengthbyLEFM_i = min(lengthbyLEFM_i);
    lengthbyLEFM = cat(1,lengthbyLEFM,lengthbyLEFM_i);
    
    
    AvgPhediStruct  = phedi_averagePhedi(DataStruct);
    
    figure(CollpsFigColor); hold on;
    [~,ci] = min(abs(DataStruct.PhediData.Cf-Cf_ref_vec));
    plot(AvgPhediStruct.X_mean,AvgPhediStruct.Vel_mean_x/(2*FWHM_Y),'Color',Cf_ref_color(ci,:),...
        'DisplayName',[DataStruct.BigPicRotStruct.details,' Cf=',num2str(DataStruct.PhediData.Cf)])
    
end


for i = 1:length(relevantEvents)
    DataStruct = phedi_smoothByfreqFilter(PhediStructCell{relevantEvents(i)});
    Movie_phedi_plotWithSg(DataStruct,'01',{'avgPhedi'});
    hold on; plot([Loc2Left_vec(i) Loc2Right_vec(i)],FWHM_Y_vec(i)*[1 1],'o','LineWidth',2)
    pause(1);
end

figure(LengthFig); hold on;
plot(Cf_vec,FWHM_vec,'*','LineWidth',2,'DisplayName',['measured FWHM', DataStruct.ExperimentData.ExpHour]);
plot(Cf_vec,abs(Loc2Right_vec),'*','LineWidth',2,'DisplayName',['right side length', DataStruct.ExperimentData.ExpHour]);
plot(Cf_vec,abs(Loc2Left_vec),'*','LineWidth',2,'DisplayName',['left side length', DataStruct.ExperimentData.ExpHour]);
plot(Cf_vec,abs(lengthbyLEFM),'*','LineWidth',2,'DisplayName',['length by LEFM cutoff', DataStruct.ExperimentData.ExpHour]);
plot([Cr Cr],[0 0.2],'LineWidth',1.5,'DisplayName','C_R');
title(['lengths ',DataStruct.ExperimentData.ExpDate,' ',DataStruct.ExperimentData.ExpHour]);

figure(CfFig); hold on;
plot(Cf_vec,2*FWHM_Y_vec,'*');
% plot([Cr Cr],[0 1],'LineWidth',1.5,'DisplayName','C_R');
title(['peak velocity ',DataStruct.ExperimentData.ExpDate,' ',DataStruct.ExperimentData.ExpHour]);



