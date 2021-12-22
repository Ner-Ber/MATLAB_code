% clc
% close all

%-----parameters

frame_crop = event_frame-500:event_frame+1500;
spatial_crop=600:900;
row=1;                                  % choose row to analize; row 9 is the mean of the 8 rows
pointsNum=1;            % choose number of points around max correlation in order to match a parabola, for a sub-pixel resolution


RowOverTime_forCorr=RowOverTimeCell{row}(frame_crop,spatial_crop);

% RowOverTime_forCorr = -RowOverTime_forCorr;

%-----Normalize each row before correlation check

% RowOverTimeShifted = RowOverTime_forCorr - repmat(min(RowOverTime_forCorr,[],2),1,size(RowOverTime_forCorr,2));
% RowOverTimeNormed_forCorr = RowOverTimeShifted./repmat(max(RowOverTimeShifted,[],2),1,size(RowOverTime_forCorr,2));
% RowOverTimeNormed_forCorr(isnan(RowOverTimeNormed_forCorr)) = 0;
RowOverTimeNormed_forCorr = RowOverTime_forCorr;                    % check correlations without normalizing



%-----correlation check

corr_cell={};
max_corr=[];
max_of_parabola_location=[];

figure;hold all
for i=1:size(RowOverTimeNormed_forCorr,1)
    
	%----correlation between each frame with the first one
    [r,lags] = xcorr(RowOverTimeNormed_forCorr(1,:),RowOverTimeNormed_forCorr(i,:));
    [f,xx,yy]=fitCorr(r',lags',pointsNum);
    
    corr_cell{i}=[r' lags'];
    max_of_parabola_location(i)=xx(find(yy==max(yy),1));
    max_corr(i) = find(r==max(r),1);
    
    
    
    %----correlation between each frame with the one before
%     if (i==1)
%         max_corr(1)=0;
%         max_of_parabola_location(1)=0;
%     else
% %         corrSamp = RowOverTimeNormed_forCorr;
% %         corrSamp = corrSamp(i,:);
% %         corrSamp([1:410,465:end]) = 0;
% %         [r,lags] = xcorr(RowOverTimeNormed_forCorr(i-1,:),corrSamp);
%         [r,lags] = xcorr(RowOverTimeNormed_forCorr(i-1,:),RowOverTimeNormed_forCorr(i,:));
%         [f,xx,yy]=fitCorr(r',lags',pointsNum);
%         
%         corr_cell{i}=[r' lags'];
%         max_of_parabola_location(i)=max_of_parabola_location(i-1)+xx(find(yy==max(yy),1));
%         [~, I] = max(r);
%         max_corr(i) = max_corr(i-1)+r(I);
%         
%     end
    
    
    %--- plot correlation functions
    if ~mod(i,50)
        %         plot(lags,r./max(r),'-','DisplayName',num2str(i))
        plot(lags,r,'-','DisplayName',num2str(i))
    end
end

max_corr = max_corr-max_corr(1);

%% graphs

shift=find(RowOverTime_forCorr(1,:)==max(RowOverTime_forCorr(1,:)));
Intensity_traj=-max_corr+shift(1);

figure;imagesc(RowOverTime_forCorr);title(['event ' num2str(eventNum) '- Row ' num2str(row)])
hold on;plot(Intensity_traj,1:size(RowOverTime_forCorr,1))


figure;hold all
plot(time(frame_crop),max_corr.*pixelsize,'-','DisplayName','MaxCorr')
plot(time(frame_crop),max_of_parabola_location.*pixelsize,'-','DisplayName','MaxofParabola')
xlim([-2 4])
xlabel('Time [msec]')
ylabel('Max Correlation [\mum]')
title('MaxCorr position')


part_vel=diff(max_of_parabola_location).*pixelsize.*IDT_Meta.FrameRate;
time_vel=time(frame_crop(1:end-1))+0.5*1/IDT_Meta.FrameRate*1e3;

figure;hold all
plot(time_vel,part_vel,'-')
plot(time_vel,smooth(part_vel),'-')
xlim([-2 4])
xlabel('Time [msec]')
ylabel('Velocity [\mum/s]')
title('MaxCorr velocity')



% Intensity=[];
Intensity_sum=[];
for i = 1:length(max_corr)
%     Intensity(i)=RowOverTime_forCorr(i,Intensity_traj(i));
    Intensity_sum(i)=sum(RowOverTime_forCorr(i,:));
end

figure;hold all
% plot(time(frame_crop),Intensity./max(Intensity))
plot(time(frame_crop),Intensity_sum./max(Intensity_sum))
xlim([-2 4])
xlabel('Time [msec]')
ylabel('Intensity')
title('Intensity of Row')