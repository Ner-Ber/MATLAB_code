function Movie_phedi_plotWithSg_23102018(DataStruct,varargin)
% Movie_phedi_plotWithSg(DataStruct,plotDamaged,PlotSeparate, plotNames, reduceSin)
%
%Movie_phedi_plotWithSg will plot phedi tracking given the output of the
%function 'Movie_phedi_from_folder_2_data';
%
% VARARGIN = DAFAULTS:
% plotDamaged=1                                         (yes, plot them)
% PlotSeparate='sub'                                    (for seperate plots enter 'sep')
% plotNames={'ROT','loc','vel','phase','locSpace','velSpace','phediUxx'}
%                                                   This is a cell array
%                                                   containing names
%                                                   of plot you want.
%                                                   relevant only in 'sep' mode.

%% set defaults
[plotDamaged,PlotSeparate, plotNames, reduceSin] = setDefaults4function(varargin,...
    1,'sub', {'ROT','loc','vel','phase','locSpace','velSpace','phediUxx'},1);

if strcmp(PlotSeparate,'sep')
    PlotSeparate=1;
else
    PlotSeparate=0;
end

%% general stuff
AsperityData = DataStruct.AsperityData;
SgData = DataStruct.SgData;
PhediData = DataStruct.PhediData;

scratchRegionMeters = [0.07 0.09];
relevantSGs = find(scratchRegionMeters(1)<SgData.x_sg*1e-3 & SgData.x_sg*1e-3<scratchRegionMeters(2));
SgColors = MyVaryColor(size(SgData.Uxx,2));
if isfield(PhediData,'Cf')
    Cf = PhediData.Cf;
else
    error('no Cf in data');
end


%     PhediLocation = -analyzePhediCell.PhediLocation;
PhediLocation = PhediData.PhediLocation;
Event = DataStruct.ExperimentData.eventNum;
FigColors = MyVaryColor(size(PhediLocation,2));
if ~PlotSeparate
    figure('Name',['Event ',num2str(Event)]);
end
LGND = cellfun(@num2str,num2cell(PhediData.measuredPhedisFromPlot),'UniformOutput',0);


center_of_mass = DataStruct.AsperityData.center_of_mass;
res = DataStruct.ExperimentData.res;
CM_time = DataStruct.AsperityData.timeCount;

%% plot the heat map with phedis
if ~(PlotSeparate&&(~any(strcmp(plotNames,'ROT'))))
    if PlotSeparate
        figure('Name',['Event ',num2str(Event),' RowOver w/ Phedis']);
    else
        subplot(2,2,1);
    end
    hold on;
    shiftByLocation = PhediData.PhotoLocation-mean(AsperityData.spatialVec);
    xImagesc = [AsperityData.spatialVec(1) AsperityData.spatialVec(end)]+shiftByLocation;
    yImagesc = [AsperityData.timeCount(1), AsperityData.timeCount(end)];
    imagesc(xImagesc,yImagesc,AsperityData.RowOverTime);
    for i = 1:size(PhediLocation,2)
        if PhediData.slopeIncline(i)>0
            plot(PhediLocation(:,i)+PhediData.measuredPhedisFromPlot(i)+shiftByLocation,...
                PhediData.timeVec,...
                '*-','Color',FigColors(i,:));
        elseif PhediData.slopeIncline(i)<0
            plot(PhediLocation(:,i)+PhediData.measuredPhedisFromPlot(i)+shiftByLocation,...
                PhediData.timeVec,...
                '.-','Color',FigColors(i,:));
        else
            if plotDamaged
                plot(PhediLocation(:,i)+PhediData.measuredPhedisFromPlot(i)+shiftByLocation,...
                    PhediData.timeVec,...
                    'r--');
            end
        end
    end
    ylabel('tims [s]');
    xlabel('approx. location on block [m]');
    title({'row over time',[DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]});
    xlim(xImagesc); ylim(yImagesc);
    %         LGND = cellfun(@num2str,num2cell(analyzePhediCell.measuredPhedisFromPlot),'UniformOutput',0);
    legend(LGND,'Location','eastoutside');
    legend('off');
    
    %--- add center of mass trajectories
    plot(center_of_mass/res+shiftByLocation,CM_time);
    
end
%% plot phedi location with sg
Phedis2Zero = bsxfun(@minus,PhediLocation,mean(PhediLocation(1:10,:),'omitnan'));
if ~(PlotSeparate&&(~any(strcmp(plotNames,'loc'))))
    if PlotSeparate
        figure('Name',['Event ',num2str(Event),' Phedis loc w/ Uxx']);
    else
        subplot(2,2,2);
    end
    
    hold on;
    %--- t-t_tip:
%     for i = 1:size(PhediLocation,2)
%         if PhediData.slopeIncline(i)>0
%             plot(PhediData.t_mins_t_tip(:,i),Phedis2Zero(:,i),...
%                 '*-','Color',FigColors(i,:));
%         elseif PhediData.slopeIncline(:,i)<0
%             plot(PhediData.t_mins_t_tip(i),Phedis2Zero(:,i),...
%                 '.-','Color',FigColors(i,:));
%         else
%             if plotDamaged
%                 plot(PhediData.t_mins_t_tip(:,i),Phedis2Zero(:,i),...
%                     'r--');
%             end
%         end
%     end
    
    %--- original time vector:
    for i = 1:size(PhediLocation,2)
        if PhediData.slopeIncline(i)>0
            plot(PhediData.timeVec,Phedis2Zero(:,i),...
                '*-','Color',FigColors(i,:));
        elseif PhediData.slopeIncline(i)<0
            plot(PhediData.timeVec,Phedis2Zero(:,i),...
                '.-','Color',FigColors(i,:));
        else
            if plotDamaged
                plot(PhediData.timeVec,Phedis2Zero(:,i),...
                    'r--');
            end
        end
    end
    
    
    %--- add cm location
%     center_of_mass_zero = bsxfun(@minus,center_of_mass,mean(center_of_mass(10:100,:),1));
%     hold on;
%     plot(CM_time,center_of_mass_zero/res,'LineWidth',1.5);

    xlabel('tims [s]');
    ylabel('location of phedi [m]');
    title('location of phedis');
    
    %--- add sg
    yyaxis right
    for i = relevantSGs
        plot(SgData.t/1000,...
            SgData.Uxx(:,i)-SgData.Uxx(1,i),...
            '.-','Color',SgColors(i,:));
    end
    ylabel('Uxx');
    LGND = [LGND;cellfun(@num2str,num2cell(SgData.x_sg(relevantSGs)'),'UniformOutput',0)];
    legend(LGND,'Location','eastoutside');
    legend('off');
    title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]);
    xlim([PhediData.timeVec(1), PhediData.timeVec(end)]);
    hold off;
end


%% plot phedi location as function of location
Phedis2Zero = bsxfun(@minus,PhediLocation,mean(PhediLocation(1:10,:),'omitnan'));
if ~(PlotSeparate&&(~any(strcmp(plotNames,'locSpace'))))
    figure('Name',['Event ',num2str(Event),' Phedis loc Vs. location']);
    
    %--- plot with reduced Sin
    if reduceSin
        PhediLocReduced = phedi_reduceSinFromLocation(PhediData);
    else
        PhediLocReduced = Phedis2Zero;
    end
    
    if isfield(PhediData,'t_mins_t_tip')
        t_mins_t_tip = PhediData.t_mins_t_tip;
    else
        t_mins_t_tip = bsxfun(@minus,PhediData.timeVec,PhediData.t_tips(:)');
    end
    
    hold on;
    %--- using differential velocity:
%     for i = 1:size(PhediLocation,2)
%         if PhediData.slopeIncline(i)>0
%             plot(PhediData.x_mins_x_tip(:,i),PhediLocReduced(:,i),...
%                 '*-','Color',FigColors(i,:));
%         elseif PhediData.slopeIncline(i)<0
%             plot(PhediData.x_mins_x_tip(:,i),PhediLocReduced(:,i),...
%                 '.-','Color',FigColors(i,:));
%         else
%             if plotDamaged
%                 plot(PhediData.x_mins_x_tip(:,i),PhediLocReduced(:,i),...
%                     'r--');
%             end
%         end
%     end
    
    %--- using constant velocity
    for i = 1:size(PhediLocation,2)
        if PhediData.slopeIncline(i)>0
            plot(-Cf*t_mins_t_tip(:,i),PhediLocReduced(:,i),...
                '*-','Color',FigColors(i,:));
        elseif PhediData.slopeIncline(i)<0
            plot(-Cf*t_mins_t_tip(:,i),PhediLocReduced(:,i),...
                '.-','Color',FigColors(i,:));
        else
            if plotDamaged
                plot(-Cf*t_mins_t_tip(:,i),PhediLocReduced(:,i),...
                    'r--');
            end
        end
    end

    xlabel('x-x_{tip} [m]');
    ylabel('location of phedi [m]');
    title('location of phedis');
    
    LGND = cellfun(@num2str,num2cell(PhediData.measuredPhedisFromPlot),'UniformOutput',0);
    legend(LGND,'Location','eastoutside');
    legend('off');
    title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]);
    xlim([-0.05 0.05]);
    hold off;
end
%% plot speed from phedis and sg
if ~(PlotSeparate&&(~any(strcmp(plotNames,'vel'))))
    %--- plot speed of phedis:
    if PlotSeparate
        figure('Name',['Event ',num2str(Event),' velocity']);
    else
        subplot(2,2,3);
    end
    hold on;
    
    originlTime=1;
    if originlTime
        [PhediVelocity, timeVec_4vel] = deal(PhediData.PhediVelocity,PhediData.timeVec_4vel);
        for i = 1:size(Phedis2Zero,2)
            if PhediData.slopeIncline(i)>0
                plot(timeVec_4vel,PhediVelocity(:,i),...
                    '*-','Color',FigColors(i,:));
            elseif PhediData.slopeIncline(i)<0
                plot(timeVec_4vel,PhediVelocity(:,i),...
                    '.-','Color',FigColors(i,:));
            else
                if plotDamaged
                    plot(timeVec_4vel,PhediVelocity(:,i),...
                        'r--');
                end
            end
        end
    else
        [PhediVelocity, timeVec_4vel] = deal(PhediData.PhediVelocity,PhediData.t_mins_t_tip_4vel);
        for i = 1:size(Phedis2Zero,2)
            if PhediData.slopeIncline(i)>0
                plot(timeVec_4vel(:,i),PhediVelocity(:,i),...
                    '*-','Color',FigColors(i,:));
            elseif PhediData.slopeIncline(i)<0
                plot(timeVec_4vel(:,i),PhediVelocity(:,i),...
                    '.-','Color',FigColors(i,:));
            else
                if plotDamaged
                    plot(timeVec_4vel(:,i),PhediVelocity(:,i),...
                        'r--');
                end
            end
        end
    end
    
    

    legendLogical = PhediData.slopeIncline<0;
    xlabel('t-t_{tip} [s]');
    ylabel('velocity of phedi [m/s]');
    ylim([-2, 2]);
    
    %--- plot particle velocity from sg
    for i = relevantSGs
        %--- particle velocity from sg:
        CurrentParticleVelocity = -SgData.v_sg(i)*(SgData.Uxx(:,i)-mean(SgData.Uxx(1:100,i)))*1e-3;
        %--- find t_tip:
        if originlTime
            SG_timeVec = (SgData.t./1000);
        else
            [~,IdxMin] = min(abs(DataStruct.BigPicRotStruct.x-SgData.x_sg(i)/1000));
            SG_t_tip = DataStruct.BigPicRotStruct.frontTime_interp(IdxMin)/DataStruct.BigPicRotStruct.fps;
            SG_timeVec = (SgData.t./1000)-SG_t_tip;
        end
        plot(SG_timeVec,CurrentParticleVelocity,...
            '.-','Color',rgb('DarkGreen'),'MarkerSize',10,'LineWidth',1);
%             '.-','Color',SgColors(i,:));
    end
    ylabel('particle velocity (m/s)');
    legend(LGND,'Location','eastoutside');
    legend('off');
    xlim([-1 1]*2e-4);
    hold off;
    title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]);
end
%% plot speed as function of location (phase space plot)
if ~(PlotSeparate&&(~any(strcmp(plotNames,'phase'))))
    if PlotSeparate
        figure('Name',['Event ',num2str(Event),' phase space']);
    else
        subplot(2,2,4);
    end
    hold on;
    [PhediVelocity, ~] = phedi_calcVelocity(Phedis2Zero,PhediData.timeVec);
    for i = 1:size(Phedis2Zero,2)
        if PhediData.slopeIncline(i)>0
            PhediLocationInterpolated = movmean(Phedis2Zero(:,i),2,'Endpoints','discard');
            
            plot(PhediLocationInterpolated,PhediVelocity(:,i),...
                '*-','Color',FigColors(i,:));
        elseif PhediData.slopeIncline(i)<0
            PhediLocationInterpolated = movmean(Phedis2Zero(:,i),2,'Endpoints','discard');
            
            plot(PhediLocationInterpolated,PhediVelocity(:,i),...
                '.-','Color',FigColors(i,:));
        else
            if plotDamaged
                PhediLocationInterpolated = movmean(Phedis2Zero(:,i),2,'Endpoints','discard');
                
                plot(PhediLocationInterpolated,PhediVelocity(:,i),...
                    'r.-');
            end
        end
    end
    legendLogical = PhediData.slopeIncline<0;
    legend(LGND,'Location','eastoutside');
    legend('off');        xlabel('relative location [m]');
    ylabel('velocity of phedi [m/s]');
%     title('location=f(velocity)');
    title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]);
    
end

%% plot particle velocity is space with cohesive zone model
if ~(PlotSeparate&&(~any(strcmp(plotNames,'velSpace'))))
    figure('Name',['Event ',num2str(Event),' vel in space+ cohesive']);
    hold on;
    [PhediVelocity, spaceVec_4vel] = deal(PhediData.PhediVelocity,PhediData.x_mins_x_tip_4vel);
    
    for i = 1:size(Phedis2Zero,2)
        if PhediData.slopeIncline(i)>0
            plot(spaceVec_4vel(:,i),PhediVelocity(:,i),...
                '*-','Color',FigColors(i,:));
        elseif PhediData.slopeIncline(i)<0
            plot(spaceVec_4vel(:,i),PhediVelocity(:,i),...
                '.-','Color',FigColors(i,:));
        else
            if plotDamaged
                plot(spaceVec_4vel(:,i),PhediVelocity(:,i),...
                    'r--');
            end
        end
    end
    
    
    %--- plot particle velocity from sg
    for i = relevantSGs
        %--- particle velocity from sg:
        CurrentParticleVelocity = -SgData.v_sg(i)*(SgData.Uxx(:,i)-mean(SgData.Uxx(1:100,i)))*1e-3;
        %--- find t_tip:
        [~,IdxMin] = min(abs(DataStruct.BigPicRotStruct.x-SgData.x_sg(i)/1000));
        SG_t_tip = DataStruct.BigPicRotStruct.frontTime_interp(IdxMin)/DataStruct.BigPicRotStruct.fps;
        SG_timeVec = (SgData.t./1000)-SG_t_tip;
        SG_x_mins_x_tip = -SG_timeVec*SgData.v_sg(i);
        plot(SG_x_mins_x_tip,CurrentParticleVelocity,...
            '.-','Color',SgColors(i,:));
    end
    
    %--- plot particle velocity at interface from cohesive zone model
    if isfield(DataStruct,'CohesiveModelStruct')
        sol=DataStruct.CohesiveModelStruct;
        x_4_Ux_dot = movmean(sol.x,2,'Endpoints','discard');
        t = -sol.x/Cf;
        Ux_dot = diff(sol.Ux)./diff(t);
        plot(x_4_Ux_dot,Ux_dot,'Color',[1 0 1],'LineWidth',2);
    else
        plot(nan,nan,'Color',[1 0 1],'LineWidth',2);
    end
    %--- make plot pretty
    xlabel('x-x_{tip} [m]');
    ylabel('velocity [m/s]');
    legend(cat(1,LGND,'cohesive model at interface'),'Location','eastoutside');
    legend('off');
    title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)])
    hold off;
    
    
end

%% plot strain calculated from phedis
if ~(PlotSeparate&&(~any(strcmp(plotNames,'phediUxx'))))
    phedisLoc4Uxx = bsxfun(@plus, PhediData.PhediLocation,PhediData.measuredPhedisFromPlot');
    % phedisNums = [4 6 8 10 12 14 16 18 22 24 26 30];
%     phedisSlopesFiction = zeros(1,size(PhediData.PhediLocation,2));
    % phedisSlopesFiction(phedisNums) = 1;
%     phedisSlopesFiction(1:length(phedisSlopesFiction)) = 1;
    
%     [UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx(phedisLoc4Uxx,[],[],0.15e-3,phedisSlopesFiction,1);
    [UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx_ofPairs(phedisLoc4Uxx);
    
    locations4Legend = cellfun(@num2str,num2cell(round(initial_locations,5)),'UniformOutput',0);
    
    figure;
    pairsLGND = {};
    hold on;
    for ii = 1:size(phediPairs,2)
        plot(PhediData.timeVec,UxxFromPhedi(:,ii),'Linewidth',1.5);
        pairsLGND{ii} = [locations4Legend{1,ii},' & ' locations4Legend{2,ii}];
    end
    UxxFromPhediMEAN = mean(UxxFromPhedi,2);
    plot(PhediData.timeVec,UxxFromPhediMEAN,'k','Linewidth',1.5);
    pairsLGND{ii+1} = 'mean Uxx';
    legend(pairsLGND,'Location','WestOutside');
    title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)])
    hold off;
end

%% main title
if ~PlotSeparate
    suptitle({['Event ',num2str(Event)],...
        ['Cf = ',num2str(Cf),' [m/s]']});
end


end