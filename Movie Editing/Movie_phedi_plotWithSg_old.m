function Movie_phedi_plotWithSg(analyzeAsperityCell, sg_data_cell, analyzePhediCell,varargin)
% Movie_phedi_plotWithSg(analyzeAsperityCell, sg_data_cell, analyzePhediCell,EventsVec,plotDamaged,PlotSeparate)
%Movie_phedi_plotWithSg will plot phedi tracking given the variables
% analyzeAsperityCell, sg_data_cell, CAM_meta_cell and analyzePhediCell in
% the workspace.
%
% defaults:
% EventsVec=1:length(analyzePhediCell)
% plotDamaged=1
% PlotSeparate='sub' (for seperate plots enter 'sep')
% plotNames={'ROT','loc','vel','phase','phediUxx'}  This is a cell array
%       containing names of plot you want. relevant only in 'sep' mode.

[EventsVec,plotDamaged,PlotSeparate, plotNames] = setDefaults4function(varargin,1:length(analyzePhediCell),1,'sub', {'ROT','loc','vel','phase','phediUxx'});
if strcmp(PlotSeparate,'sep')
    PlotSeparate=1;
else
    PlotSeparate=0;
end

for Event = EventsVec
    %--- pass to next iteration of there is no data:
    if isempty(analyzePhediCell{Event})
        continue
    end
    
    
    %     PhediLocation = -analyzePhediCell{Event}.PhediLocation;
    PhediLocation = analyzePhediCell{Event}.PhediLocation;
    FigColors = MyVaryColor(size(PhediLocation,2));
    if ~PlotSeparate
        figure('Name',['Event ',num2str(Event)]);
    end
    LGND = cellfun(@num2str,num2cell(analyzePhediCell{Event}.measuredPhedisFromPlot),'UniformOutput',0);
    
    %% plot the heat map with phedis
    if ~(PlotSeparate&&(~any(strcmp(plotNames,'ROT'))))
        if PlotSeparate
            figure('Name',['Event ',num2str(Event),' RowOver w/ Phedis']);
        else
            subplot(2,2,1);
        end
        hold on;
        %         xImagesc = -[analyzeAsperityCell{Event}.spatialVec(end) analyzeAsperityCell{Event}.spatialVec(1)];
        xImagesc = [analyzeAsperityCell{Event}.spatialVec(1) analyzeAsperityCell{Event}.spatialVec(end)];
        yImagesc = [analyzeAsperityCell{Event}.timeCount(1), analyzeAsperityCell{Event}.timeCount(end)];
        imagesc(xImagesc,yImagesc,analyzeAsperityCell{Event}.RowOverTime);
        for i = 1:size(PhediLocation,2)
            if analyzePhediCell{Event}.slopeIncline(i)>0
                plot(PhediLocation(:,i)+analyzePhediCell{Event}.measuredPhedisFromPlot(i),...
                    analyzePhediCell{Event}.timeVec,...
                    '*-','Color',FigColors(i,:));
            elseif analyzePhediCell{Event}.slopeIncline(i)<0
                plot(PhediLocation(:,i)+analyzePhediCell{Event}.measuredPhedisFromPlot(i),...
                    analyzePhediCell{Event}.timeVec,...
                    '.-','Color',FigColors(i,:));
            else
                if plotDamaged
                    plot(PhediLocation(:,i)+analyzePhediCell{Event}.measuredPhedisFromPlot(i),...
                        analyzePhediCell{Event}.timeVec,...
                        'r--');
                end
            end
        end
        ylabel('tims [s]');
        xlabel('distance [m]');
        title('row over time');
        xlim(xImagesc); ylim(yImagesc);
        %         LGND = cellfun(@num2str,num2cell(analyzePhediCell{Event}.measuredPhedisFromPlot),'UniformOutput',0);
        legend(LGND,'Location','eastoutside');
        legend('off');
    end
    %% plot phedi location with sg
    %         Phedis2Zero = PhediLocation-repmat(mean(PhediLocation(1:100,i),1),size(PhediLocation,1),size(PhediLocation,2));
    Phedis2Zero = bsxfun(@minus,PhediLocation,PhediLocation(1,:));
    
    if ~(PlotSeparate&&(~any(strcmp(plotNames,'loc'))))
        if PlotSeparate
            figure('Name',['Event ',num2str(Event),' Phedis loc w/ Uxx']);
        else
            subplot(2,2,2);
        end
        
        hold on;
        for i = 1:size(PhediLocation,2)
            if analyzePhediCell{Event}.slopeIncline(i)>0
                plot(analyzePhediCell{Event}.timeVec,Phedis2Zero(:,i),...
                    '*-','Color',FigColors(i,:));
            elseif analyzePhediCell{Event}.slopeIncline(i)<0
                plot(analyzePhediCell{Event}.timeVec,Phedis2Zero(:,i),...
                    '.-','Color',FigColors(i,:));
            else
                if plotDamaged
                    plot(analyzePhediCell{Event}.timeVec,Phedis2Zero(:,i),...
                        'r--');
                end
            end
        end
        xlabel('tims [s]');
        ylabel('location of phedi [m]');
        title('location of phedis');
        
        %--- add sg
        scratchRegionMeters = [0.07 0.09];
        relevantSGs = find(scratchRegionMeters(1)<sg_data_cell{Event}.x_sg*1e-3 & sg_data_cell{Event}.x_sg*1e-3<scratchRegionMeters(2));
        SgColors = MyVaryColor(size(sg_data_cell{Event}.Uxx,2));
        yyaxis right
        for i = relevantSGs
            plot(sg_data_cell{Event}.t/1000,...
                sg_data_cell{Event}.Uxx(:,i)-sg_data_cell{Event}.Uxx(1,i),...
                '.-','Color',SgColors(i,:));
        end
        ylabel('Uxx');
        LGND = [LGND;cellfun(@num2str,num2cell(sg_data_cell{Event}.x_sg(relevantSGs)'),'UniformOutput',0)];
        legend(LGND,'Location','eastoutside');
        legend('off');
        if isfield(sg_data_cell{Event},'Cf')
            Cf = sg_data_cell{Event}.Cf;
        else
            Cf = sg_calc_velocity_from_sg(sg_data_cell{Event}, relevantSGs);
        end
        title('phedi location with sg');
        xlim([-1 1]*2e-4);
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
        
        [PhediVelocity, timeVec_4vel] = phedi_calcVelocity(Phedis2Zero,analyzePhediCell{Event}.timeVec);
        
        for i = 1:size(Phedis2Zero,2)
            if analyzePhediCell{Event}.slopeIncline(i)>0
                plot(timeVec_4vel,PhediVelocity(:,i),...
                    '*-','Color',FigColors(i,:));
            elseif analyzePhediCell{Event}.slopeIncline(i)<0
                plot(timeVec_4vel,PhediVelocity(:,i),...
                    '.-','Color',FigColors(i,:));
            else
                if plotDamaged
                    plot(timeVec_4vel,PhediVelocity(:,i),...
                        'r--');
                end
            end
        end
        
        legendLogical = analyzePhediCell{Event}.slopeIncline<0;
        xlabel('tims [s]');
        ylabel('velocity of phedi [m/s]');
        ylim([-2, 2]);
        
        %--- plot particle velocity from sg
        for i = relevantSGs
            CurrentParticleVelocity = -Cf*(sg_data_cell{Event}.Uxx(:,i)-mean(sg_data_cell{Event}.Uxx(1:100,i)))*1e-3;
            plot(sg_data_cell{Event}.t./1000,CurrentParticleVelocity,...
                '.-','Color',SgColors(i,:));
        end
        ylabel('particle velocity (m/s)');
        legend(LGND,'Location','eastoutside');
        legend('off');
        xlim([-1 1]*2e-4);
        hold off;
    end
    %% plot speed as function of location (phas espace plot)
    if ~(PlotSeparate&&(~any(strcmp(plotNames,'phase'))))
        if PlotSeparate
            figure('Name',['Event ',num2str(Event),' phase space']);
        else
            subplot(2,2,4);
        end
        hold on;
        [PhediVelocity, ~] = phedi_calcVelocity(Phedis2Zero,analyzePhediCell{Event}.timeVec);
        for i = 1:size(Phedis2Zero,2)
            if analyzePhediCell{Event}.slopeIncline(i)>0
                PhediLocationInterpolated = movmean(Phedis2Zero(:,i),2,'Endpoints','discard');
                
                plot(PhediLocationInterpolated,PhediVelocity(:,i),...
                    '*-','Color',FigColors(i,:));
            elseif analyzePhediCell{Event}.slopeIncline(i)<0
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
        legendLogical = analyzePhediCell{Event}.slopeIncline<0;
        legend(LGND,'Location','eastoutside');
        legend('off');        xlabel('relative location [m]');
        ylabel('velocity of phedi [m/s]');
        title('location=f(velocity)');
        
    end
    
    %% plot strain calculated from phedis
    if ~(PlotSeparate&&(~any(strcmp(plotNames,'phediUxx'))))
        phedisLoc4Uxx = bsxfun(@plus, analyzePhediCell{Event}.PhediLocation,analyzePhediCell{Event}.measuredPhedisFromPlot');
        % phedisNums = [4 6 8 10 12 14 16 18 22 24 26 30];
        phedisSlopesFiction = zeros(1,size(analyzePhediCell{Event}.PhediLocation,2));
        % phedisSlopesFiction(phedisNums) = 1;
        phedisSlopesFiction(1:length(phedisSlopesFiction)) = 1;
        
        [UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx(phedisLoc4Uxx,[],[],0.3e-3,phedisSlopesFiction,1);
        
        locations4Legend = cellfun(@num2str,num2cell(round(initial_locations,5)),'UniformOutput',0);
        
        figure;
        pairsLGND = {};
        hold on;
        for ii = 1:size(phediPairs,2)
            plot(analyzePhediCell{Event}.timeVec,UxxFromPhedi(:,ii),'Linewidth',1.5);
            pairsLGND{ii} = [locations4Legend{1,ii},' & ' locations4Legend{2,ii}];
        end
        UxxFromPhediMEAN = mean(UxxFromPhedi,2);
        plot(analyzePhediCell{Event}.timeVec,UxxFromPhediMEAN,'k','Linewidth',1.5);
        pairsLGND{ii+1} = 'mean Uxx';
        legend(pairsLGND,'Location','WestOutside');
        hold off;
    end
    
    
    %% main title
    if ~PlotSeparate
        suptitle({['Event ',num2str(Event)],...
            ['Cf = ',num2str(Cf),' [m/s]']});
    end
end



end