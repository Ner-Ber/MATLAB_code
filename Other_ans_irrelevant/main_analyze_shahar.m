
% clear
% close all
clc

exper=my_dir;
experNum = 3;
eventNum=6;  %% choose eventNo
row=9;

folow_asperities=0;
plot_asperities_movment=0;

%% Create the Images matrix, [rows * spatial * time] 

pre_frames=1000;
post_frames=4000;
% [Images,Cam_Meta] = IDT_ReadImages(exper{experNum},eventNum,pre_frames,1,post_frames,0);   % IDT
[Images,Cam_Meta] = photronReadImages(exper{experNum},eventNum,pre_frames,1,post_frames,0);   %Photron

%% Make The Row-Over-Time Pictures, each cell for a different row

Rows=1:size(Images,1);                 %% choose rows to look at
Mean_RowOverTime=mean(Images,1);
Mean_RowOverTime=permute(Mean_RowOverTime,[3 2 1]);

RowOverTimeCell = {};
for rowIdx = 1:length(Rows)
    RowOverTime = Images(rowIdx,:,:);
    RowOverTime = permute(RowOverTime,[3 2 1]);
    RowOverTimeCell{rowIdx} = RowOverTime;
end

RowOverTimeCell{length(Rows)+1}=Mean_RowOverTime;
clear('Mean_RowOverTime');clear('RowOverTime');

%% Define the Asperities to follow

if folow_asperities==1
    
    % factor=1.3;            % factor for threshold, asperities will be choosed only if they are higher than it
    % min_gap=5;             % minimum gap between asperities
    % asperity_min_size=12;  % minimum width of asperities
    %
    %
    % [vallyPairsCell,asper_num]=define_asperities(RowOverTimeCell,Rows,factor,min_gap,asperity_min_size);
    % disp('num of asperities in each row:');
    % disp(asper_num);
    
    vallyPairsCell = my_vallyPairsCell;
end

%% Make the images ready for analizing 

if folow_asperities==1
    
    % event_frame_smallcrop=pre_frames_bigcrop+1;       %% another crop for further investigation
    % pre_frames_smallcrop=50;
    % post_frames_smallcrop=50;
    % RelevantFrames=event_frame_smallcrop-pre_frames_smallcrop:...
    %     event_frame_smallcrop+post_frames_smallcrop;
    
    % [~,RowOverTimeCropedCell, ~,...                   % just the small-crop cell
    %     ~, ~, ~] =...
    %     Movie_follow_maxima_in_row(Images,vallyPairsCell, Rows, RelevantFrames, [], 1, []);
    
    [RowOverTimeNormCell,~, ~,MarkedMatCell, AsperityDefinitionTrajCell] =...
        Movie_follow_maxima_in_row(Images,vallyPairsCell, Rows, 'all', [], 1, []);   % all other cells
end

%% Analyze the picture to get the trajectories of desired data:

if folow_asperities==1
    [maxTrajectories, CM_Trajectories, skewnessCell, varCell,center_of_massCell] =...
        Movie_analyze_asperities(MarkedMatCell, RowOverTimeCell);
end

%% make time vector

time=1:size(Images,3);
time=time/Cam_Meta.FrameRate;                     % [sec]
event_frame=pre_frames+1;
time=time-time(event_frame);

pixelsize=3e-6;      % [m]

%% Get data from SG:
sg_data=acq132_event_get_data...
    (exper{experNum},eventNum,'start','end',1,'Uxx','Uyy','Uxy','x_sg','F','N');

x_laser=105;

[c_f,relevant_sg,shifted_U]=FindRuptureVelocity(sg_data,x_laser);

%% plot Uxx

figure;plot(sg_data.t,shifted_U)
legend(num2str(sg_data.x_sg'))
legend off
title({'Uxx';['event ' num2str(eventNum) '- Row ' num2str(row)];['C_f = ' num2str(c_f) ' [m/s]']})
xlabel('Time [msec]')

xlim([-2 4])

%% Graphs:
if plot_asperities_movment==1
    row=1;
    
    %--- plot 1: Row-over-Time of a single row with the trajectories on it
    CM = CM_Trajectories{row};
    
    figure;
    % subplot(2,2,1)
    colormap jet; AsperColor = MyVaryColor(length(CM)); colormap parula;
    
    imagesc(RowOverTimeCell{row});
    hold on;
    for i = 1:length(CM)
        plot(CM{i}(:,1),CM{i}(:,2),'-','Color',AsperColor(i,:));
    end
    legend(num2str(vallyPairsCell{row}(:,1)))
    legend off
    plot(AsperityDefinitionTrajCell{row}(:,2),AsperityDefinitionTrajCell{row}(:,1),'g.');
    title(['Row' num2str(row)])
    ylim([600 1800])
    
    
    %--- plot 2: CM of all asperities in one row
    
    CM = CM_Trajectories{row};
    % subplot(2,2,2);
    figure
    hold on
    for i = 1:length(CM)
        plot(time,(-CM{i}(:,1)+CM{i}(event_frame-10,1)).*pixelsize,...
            '-','Color',AsperColor(i,:));                     %convert from pixels to mu m
    end
    % line([time(event_frame) time(event_frame)],[-5 10],'Color','k','DisplayName','event');
    title(['Center of mass - Row' num2str(row)])
    legend(num2str(vallyPairsCell{row}(:,1)))
    legend off
    xlabel('Time [msec]')
    ylabel('Shifted CM [\mum]')
    xlim([-0.5 2])
    
    
    %--- plot 3: CM Intensity of all asperities in one row
    CM = CM_Trajectories{row};
    Intensity=[];
    
    figure; hold on
    for i = 1:length(CM)
        for j=1:length(CM{i})
            if isnan(CM{i}(j,1))
                Intensity(j,i)=nan;
            else
                Intensity(j,i)=RowOverTimeCell{row}(CM{i}(j,2),round(CM{i}(j,1)));
            end
        end
        m=mean(Intensity(1:event_frame-100,i));
        plot(time,Intensity(:,i)./m,...
            '-','Color',AsperColor(i,:),'DisplayName',['Asperity ' num2str(i)]);
    end
    
    title(['CM Intensity - Row' num2str(row)])
    xlabel('Time [msec]')
    ylabel('Normalized Intensity')
    xlim([-0.5 2])
end

%% mean graphs
if plot_asperities_movment==1
    
    CM={};
    CM_reduced={};
    for row=Rows
        CM{row}=[];
        k=1;
        for i=1:length(CM_Trajectories{row})
            m=mean(CM_Trajectories{row}{i}(1:event_frame-100,1));
            CM{row}(:,i) = -(CM_Trajectories{row}{i}(:,1)-m);
            if ~isnan(CM_Trajectories{row}{i}(:,1))
                CM_reduced{row}(:,k) = -(CM_Trajectories{row}{i}(:,1)-m);
                k=k+1;
            end
            
        end
        CM_reduced{row}(:,k) = mean(CM_reduced{row},2);
    end
    
    
    row2present=4;
    
    figure;hold on
    plot(time,CM{row2present}.*pixelsize)
    legend(num2str(vallyPairsCell{row2present}(:,1)))
    legend off
    title(['Center of mass - Row' num2str(row2present)])
    xlabel('Time [msec]')
    ylabel('Shifted CM [\mum]')
    xlim([-2 4])
    
    
    MeanCM=[];
    figure; hold all
    for rowIdx=1:length(CM_reduced)
        m=CM_reduced{rowIdx}(:,size(CM_reduced{rowIdx},2));
        plot(time,m.*pixelsize,'DisplayName',['Mean - Row ' num2str(rowIdx)])
        %     if ~isnan(m)
        MeanCM(:,rowIdx)=m;
        %     end
    end
    
    plot(time,mean(MeanCM,2).*pixelsize,'DisplayName','Mean of Mean')
    
    title('Center of mass - mean')
    xlabel('Time [msec]')
    ylabel('Shifted CM [\mum]')
    xlim([-2 4])
end
   
%% other graphs
% %--- plot 4: Var of all asperities in one row
% Var = varCell{row};
% 
% % subplot(2,2,3);
% figure
% hold on
% for i = 1:size(Var,2)
%     m=mean(Var(1:event_frame-100,i));
%     plot(time,Var(:,i)-m,...
%         '-','Color',AsperColor(i,:));
% end
% % line([time(event_frame) time(event_frame)],[-0.03 0.03],'Color','k','DisplayName','event');
% title(['Variance-Row' num2str(row)])
% legend(num2str(vallyPairsCell{row}(:,1)))
% legend off
% xlabel('Time [msec]')
% ylabel('Normalized Variance')
% % xlim([-0.5 2])
% 
% 
% %--- plot 5: Skewness of all asperities in one row
% Skewness = skewnessCell{row};
% 
% % subplot(2,2,4);
% figure;
% hold on
% for i = 1:size(Skewness,2)
%     m=mean(Skewness(1:event_frame-100,i));
%     plot(time,Skewness(:,i),...
%         '-','Color',AsperColor(i,:));
% end
% line([time(event_frame) time(event_frame)],[-1 1],'Color','k','DisplayName','event');
% title(['Skewness-Row' num2str(row)])
% legend(num2str(vallyPairsCell{row}(:,1)))
% legend off
% xlabel('Time [msec]')
% ylabel('Normalized Skewness')
% % xlim([-0.5 2])

% %--- plot 6: cm of all asperities in one row
% cm = center_of_massCell{row};
% 
% % subplot(2,2,4);
% figure;
% hold on
% for i = 1:size(cm,2)
%     m=mean(cm(1:event_frame-100,i));
%     plot(time,-cm(:,i)+m,...
%         '-','Color',AsperColor(i,:));
% end
% line([time(event_frame) time(event_frame)],[-1 1],'Color','k','DisplayName','event');
% title(['cm-Row' num2str(row)])
% legend(num2str(vallyPairsCell{row}(:,1)))
% legend off
% xlabel('Time [msec]')
% % ylabel('Normalized Skewness')
% % xlim([-2 4])

