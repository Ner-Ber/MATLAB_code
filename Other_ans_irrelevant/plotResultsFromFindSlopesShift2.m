% figure; hold on;
% colormap jet;
% Current_frames = 950:1200;
% FigColors = MyVaryColor(length(Current_frames));
% my_x_axis = 1:size(RowOverTimeCell{1},2);
% for i = 1:length(Current_frames)
%     plot(my_x_axis-PhediLocation(ismember(frameCount,Current_frames(i))),...
%         RowOverTimeCell{1}(Current_frames(i),:),...
%         'Color',FigColors(i,:));
%
% end

DoPlots = 0;
preFrames = 300;
postFrames = 4000;
experNum = 4;
eventNums_vec = [5]';

analyzeAsperityCell = cell(max(eventNums_vec),1);
sg_data_cell = cell(max(eventNums_vec),1);
CAM_meta_cell = cell(max(eventNums_vec),1);

analyzePhediCell = cell(max(eventNums_vec),1);

for eventNum = eventNums_vec(:)'
    entireLoopTimer = tic;
    [RowOverTimeCell, vallyPairsCell,~, ~,...
        maxTrajectories, CM_Trajectories, skewnessCell, varCell,center_of_massCell, sg_data, ~,Cam_Meta]...
        = analyzeEnlargementAndStrain(experNum, eventNum,DoPlots, preFrames, postFrames);
    %--- save asperity data
    disp('saving asperity data...');
    analyzeAsperityStruct = struct;
    [analyzeAsperityStruct.RowOverTime,...
        analyzeAsperityStruct.vallyPairs,...
        analyzeAsperityStruct.maxTrajectories,...
        analyzeAsperityStruct.CM_Trajectories,...
        analyzeAsperityStruct.skewness,...
        analyzeAsperityStruct.var,...
        analyzeAsperityStruct.center_of_mass,...
        sg_data_cell{eventNum},...
        CAM_meta_cell{eventNum}]...
        = deal(...
        RowOverTimeCell{9},...
        vallyPairsCell{9},maxTrajectories{9}, CM_Trajectories{9},...
        skewnessCell{9}, varCell{9},center_of_massCell{9}, sg_data,Cam_Meta);
    
    analyzeAsperityCell{eventNum} = analyzeAsperityStruct;
    
    %% General stuff
    
    relevanrRowOverTime = RowOverTimeCell{9};
    res = 32/(100e-6);    % pix/m
    fps = 8e5;        % frames/s
    
    timeCount = [-preFrames:postFrames]/fps;
    spatialVec = (0:(size(RowOverTimeCell{1},2)-1))/res;
    
    %% calculate phedis
    disp('analyzing phedis...');
    % phedislocationInPix = phedisDB(2);
    smoothParameter = 0.001;
    phedislocationInPix = Movie_phedi_findTrenches(relevanrRowOverTime(270,:), smoothParameter);
    %--- eliminate phedis close to edge:
    safetyDist = 20;
    phedislocationInPix = phedislocationInPix(abs(phedislocationInPix-length(relevanrRowOverTime(270,:)))>safetyDist);
    phedislocationInPix = phedislocationInPix(phedislocationInPix>safetyDist);
    phedislocationInPix = phedislocationInPix(2:(end-1));
    calShiftTimer = tic;
    [frameCount, PhediLocationPix, slopeIncline] = Movie_phedi_calculate_continues_shift(relevanrRowOverTime, 285, 420, phedislocationInPix);
    T_shiftTimer = toc(calShiftTimer);
    disp(['tie for claculatin phedi shifts = ',num2str(T_shiftTimer)]);
    
    %--- transform units:
    PhediLocation = (PhediLocationPix-1)/res;
    measuredPhedisFromPlot = (phedislocationInPix-1)/res;
    timeVec= timeCount(frameCount);
    
    %--- save Phedi data
    disp('saving phedi data...');
    analyzePhediStruct = struct;
    [analyzePhediStruct.PhediLocation,...
        analyzePhediStruct.timeVec,...
        analyzePhediStruct.measuredPhedisFromPlot,...
        analyzePhediStruct.slopeIncline]...
        = deal(...
        PhediLocation, timeVec, measuredPhedisFromPlot, slopeIncline);
    analyzePhediCell{eventNum} = analyzePhediStruct;
    
    %% plot heat map and phedis dynamics
    PlotAll = 0;
    if PlotAll
        FigColors = MyVaryColor(size(PhediLocation,2));
        figure; hold on;
        colormap default;
        imagesc([spatialVec(1) spatialVec(end)],[timeCount(1) timeCount(end)],relevanrRowOverTime);
        colormap jet
        for i = 1:size(PhediLocation,2)
            if slopeIncline(i)>0
                plot(PhediLocation(:,i)+measuredPhedisFromPlot(i),timeVec,...
                    '*-','Color',FigColors(i,:));
            elseif slopeIncline(i)<0
                plot(PhediLocation(:,i)+measuredPhedisFromPlot(i),timeVec,...
                    '.-','Color',FigColors(i,:));
            else
                plot(PhediLocation(:,i)+measuredPhedisFromPlot(i),timeVec,...
                    'r--');
            end
        end
        LGND = cellfun(@num2str,num2cell(measuredPhedisFromPlot),'UniformOutput',0);
        LGND = cat(1,LGND,{'row over time'});
        legend(LGND,'location','eastoutside');
        %% plot only phedi
        
        figure; hold on;
        colormap jet
        for i = 1:size(PhediLocation,2)
            if slopeIncline(i)>0
                        plot(timeVec,PhediLocation(:,i),...
                            '*-','Color',FigColors(i,:));
            elseif slopeIncline(i)<0
                plot(timeVec,PhediLocation(:,i),...
                    '.-','Color',FigColors(i,:));
            else
                %         plot(timeVec,PhediLocation(:,i),...
                %             'r--');
            end
        end
        LGND = cellfun(@num2str,num2cell(measuredPhedisFromPlot),'UniformOutput',0);
        legend(LGND);
        title('phedi location as function of time');
        xlabel('time [s]');
        ylabel('relative distance [m]');
        
        %% *phedi velocity*
        
        Cf = 105;        %(m/s)
        PhediVelocity = -diff(PhediLocation,[],1)./diff(repmat(timeVec',1,size(PhediLocation,2)),[],1);
        velocityTime = movmean(timeVec,2,'Endpoints','discard');
        
        figure; hold on;
        colormap jet
        for i = 1:size(PhediVelocity,2)
            if slopeIncline(i)>0
                %         plot(velocityTime,PhediVelocity(:,i),...
                %             '*-','Color',FigColors(i,:));
            elseif slopeIncline(i)<0
                plot(velocityTime,PhediVelocity(:,i),...
                    '.-','Color',FigColors(i,:));
            else
                %         plot(velocityTime,PhediVelocity(:,i),...
                %             'r--');
            end
        end
        LGND = cellfun(@num2str,num2cell(measuredPhedisFromPlot),'UniformOutput',0);
        legend(LGND);
    end
    T_iteration = toc(entireLoopTimer);
    disp(['Elapsed time for this iteration is ',num2str(T_iteration),' seconds']);
end

%% plot crossection
% colormap jet;
% Movie_PlotRowFromHeatMap(relevanrRowOverTime,850:10:(850+1750));
