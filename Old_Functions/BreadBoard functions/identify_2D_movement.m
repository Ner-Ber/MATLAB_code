
%% old - choose video:
%%% 1 for 100 mps
%%% 2 for 250 mps
%%% 3 for 15.9__4
%%% 4 for 19.9__12
%%% 5 for 27.10_20
%%% 6 for 27.10__21
%%% 7 for 27.10__18

 
% [MoviePath,CoordinateCell,EventFrame]=Choose_video(video_number);

%% pick movie
Date='\2_11\';
video_number='4';
% video_number = num2str(video_number);

Path='C:\Users\owner\Documents\Shahar_Neri\friction_videos';
Main_Folder=[Path Date video_number];

%% define starting frame and regions

Name_of_video_folder=['\Measure' video_number '_C001H001S0001'];
figures_folder_path=[Main_Folder Name_of_video_folder];
video_struct = [Main_Folder '\analyzeDataStruct.mat']; 

%# if this video was never analized before you need to define the starting...
...frame of the event and the interesting regions you want to follow
%# you do this with this fuction and the instructions in it
switch exist(video_struct, 'file')
    case 0
        [analyzeDataStruct] = AnalyzeSingleMovie(figures_folder_path);
        save(video_struct, 'analyzeDataStruct');
end

load(video_struct);
MoviePath = analyzeDataStruct.MoviePath;
EventFrame_num = analyzeDataStruct.EventFrame;
CoordinateCell = analyzeDataStruct.CoordinateCell;%% other parameters

%# when no movie_name or movieType are inserted no movie will be generated, only frames.
%# choose if to fix blinking (1) or not (0)
fix_blinking=1;
movie_name = '';
movieType = '';

%---measurment parameters
ImageData = get_cih_data(MoviePath);
dt = 1/ImageData.RecordRate_fps;
pixel_resolution = 3e-6;        %pixel capture this much meters

%---relevant time crop
pre_frames=30;
post_frames=50; %length(relevant_frames);
relevant_frames = (EventFrame_num-pre_frames):(EventFrame_num+post_frames);

%---interpolation
InterpMethod = 'cubic';
InterpFactor = 1;

%---Aspereties movies
make_asperities_movies=0;

%% find movement of displacement sensors

calibration_file=[Main_Folder '\calibration' video_number '.csv'];
new_calibration_file_name=[calibration_file '_new.mat'];

%# check if the sensors in this folder are calibrated.
%# if not - calibrate them using the best (closest in time) calibration measurment
%# before running the code you should put a calibration measurment in the folder...
...and name it 'calibration+number_of folder'
    
switch exist(new_calibration_file_name, 'file')
    case 0
        [~,~,~,~,~,new_calibration_file] = Calibration(calibration_file,'csv');
end

%# use the calibration to define the displacement sensors movement during the event
event_file=[Main_Folder '\Measure' video_number '.csv'];
[D1,D2,T,~,~,~,V_rupture] = Find_movement_of_sensors(new_calibration_file,event_file,'csv');

%% generate RGB intensity movie
[fixed_frames, fixed_frames_RGB, all_frames] =...
    generateFixedMovie(MoviePath, movie_name, movieType,fix_blinking);

%% analyze frames
[Ylength, Xlength, Tlength] = size(fixed_frames);

%---define time scale and crop in time:
frameScale = 1:Tlength;
timeScale = (frameScale-frameScale(1)+ImageData.StartFrame)./ImageData.RecordRate_fps;

frame_range = frameScale(relevant_frames);
time_range = timeScale(relevant_frames);

fixed_frames_cr = fixed_frames(:,:,relevant_frames);
fixed_frames_RGB_cr = fixed_frames_RGB(:,:,relevant_frames);
all_frames_cr = all_frames(:,:,relevant_frames);

%% find peaks after interpolation

[~, CoorOrder] = sort(cellfun(@(a) a(1), CoordinateCell));
CoordinateCell = CoordinateCell(CoorOrder);

%---create mesh
[Xgrid Ygrid] = meshgrid(1:Xlength, 1:Ylength);
% [XgridI YgridI] = meshgrid(linspace(1,Xlength,Xlength*InterpFactor),...
%     linspace(1,Ylength,Ylength*InterpFactor));


%---process each region seperatly
XmeshesCell = cell(size(CoordinateCell));                                              % cell aray to save all meshes in X coordinate
YmeshesCell = cell(size(CoordinateCell));                                              % cell aray to save all meshes in Y coordinate
interppesFramesCell = cell(size(CoordinateCell));                          % cell aray to save all interpolated movies. each cell contains a cropped movie
Extremas = cell(length(CoordinateCell),length(frame_range));   % cell array to save all data about extremas in certain frame. ...
                                                                                                                                             %each row is a different region define in coordinatesCell and each column is a time step
                                                                                                                                             
for regionIdx = 1:length(CoordinateCell)
    x_1 = min(CoordinateCell{regionIdx}([1,3]));
    x_2 = max(CoordinateCell{regionIdx}([1,3]));
%     y_1 = min(CoordinateCell{regionIdx}([2,4]));
%     y_2 = max(CoordinateCell{regionIdx}([2,4]));
    y_1=1;
    y_2=8;
    currentRegionFrames = fixed_frames_cr(y_1:y_2,x_1:x_2,:);
    
    %---interpolate frames
    currentXgrid = Xgrid(y_1:y_2,x_1:x_2);
    currentYgrid = Ygrid(y_1:y_2,x_1:x_2);
    [currentXgridI, currentYgridI] = meshgrid(linspace(x_1,x_2,(x_2-x_1)*InterpFactor),...
    linspace(y_1,y_2,(y_2-y_1)*InterpFactor));
    
    interpped_frames = zeros(size(currentXgridI,1),size(currentXgridI,2),length(frame_range));
    for timeStep = 1:length(frame_range)
        interpped_frames(:,:,timeStep) = interp2(currentXgrid, currentYgrid,...
            currentRegionFrames(:,:,timeStep), currentXgridI, currentYgridI, InterpMethod);
        %---find extremas
        [frameMAX,IMAX,~,~] = extrema2(interpped_frames(:,:,timeStep));
        Extremas{regionIdx,timeStep} = [frameMAX,IMAX];
    end
 
    %---save meshes and interpolated frames
    XmeshesCell{regionIdx} = currentXgridI;
    YmeshesCell{regionIdx} = currentYgridI;
    interppesFramesCell{regionIdx} = interpped_frames;
 
end

%% make movies by regions
switch make_asperities_movies
    case 1
        AsperitiesMovies=cell(size(CoordinateCell));
        Moviefigures=cell(size(CoordinateCell));
        for regionIdx = 1:length(CoordinateCell)
            [AsperitiesMovies{regionIdx},Moviefigures{regionIdx}] = meshMovie2(interppesFramesCell{regionIdx}...
            ,XmeshesCell{regionIdx},YmeshesCell{regionIdx});
            close
        end
end

%% follow mass center

% Makes matrixes of the mass center at each frame and each region
Xcenter = zeros(length(CoordinateCell),length(frame_range));
Ycenter = zeros(length(CoordinateCell),length(frame_range));
 for regionIdx = 1:length(CoordinateCell)
     for timeStep = 1:length(frame_range)
         [X_c,Y_c]=Center_of_Mass(interppesFramesCell{regionIdx}(:,:,timeStep),InterpFactor);
         Xcenter(regionIdx,timeStep)=X_c+CoordinateCell{regionIdx}(1)-1;
         Ycenter(regionIdx,timeStep)=Y_c+CoordinateCell{regionIdx}(2)-1;
     end
 end

%% Units

Xcenter=Xcenter*pixel_resolution;   %%pixel to meter
Ycenter=Ycenter*pixel_resolution;   
Time=[1:length(relevant_frames)]*dt;  %%frame to second  
 
%% follow peaks

% % --- find location of peaks over time
% MaxIdxs = cellfun(@(a) a(1,2), Extremas);       % creates and array containing of index of maxima in each frame. 
%                                                                                            %each row is a different region defined in coordinatesCell and each column is a time step
% MaxXcoor = zeros(size(MaxIdxs));
% MaxYcoor = zeros(size(MaxIdxs));
% for timeStep = 1:length(frame_range)        %iterate on all time steps (frames)
%     MaxXcoor(:,timeStep) = cellfun(@(a, b) a(b), XmeshesCell, num2cell(MaxIdxs(:,timeStep)));
%     MaxYcoor(:,timeStep) = cellfun(@(a, b) a(b), YmeshesCell, num2cell(MaxIdxs(:,timeStep)));
% end

%% follow speed
% Xspeed = diff(MaxXcoor, 1, 2);
% Yspeed = diff(MaxYcoor, 1, 2);
% speedMag = sqrt(Xspeed.^2 + Yspeed.^2);

%% plots we are currently using

ColorSet = MyVaryColor(length(CoordinateCell));

%---plot1 - superposition of the movement in X direction
figure;
hold all;
Xcenter_normalized=zeros(size(Xcenter));
meanVecX = [];
for regionIdx=1:length(CoordinateCell)
    mean_position=mean(Xcenter(regionIdx,1:pre_frames));
    Xcenter_normalized(regionIdx,:)=Xcenter(regionIdx,:)-mean_position;
    plot(Time,Xcenter_normalized(regionIdx,:),...
        '-','Color',ColorSet(regionIdx,:),'LineWidth',2);
    meanVecX = cat(1, meanVecX, mean_position/pixel_resolution);
end

title({'Superposition of events - Mass center on X axis',...
    ['date:   ',Date,'      video:   ',video_number]});
legend(num2str(meanVecX),'Location','NorthEastOutside');
xlabel('Time [s]');
ylabel('normalized X_c [m]' );
saveas(gcf,[Main_Folder,'\','X_movement.png']);


%---plot2 - superposition of the movement in Y direction
figure;
hold all;
Ycenter_normalized=zeros(size(Ycenter));
meanVecY = [];
for regionIdx=1:length(CoordinateCell)
    mean_position=mean(Ycenter(regionIdx,1:pre_frames));
    Ycenter_normalized(regionIdx,:)=Ycenter(regionIdx,:)-mean_position;
    plot(Time,Ycenter_normalized(regionIdx,:),...
        '-','Color',ColorSet(regionIdx,:),'LineWidth',2);
    meanVecY = cat(1, meanVecY, mean_position/pixel_resolution);
end
title({'Superposition of events - Mass center on Y axis',...
    ['date:   ',Date,'      video:   ',video_number]});
legend(num2str(meanVecX),'Location','NorthEastOutside');
xlabel('Time [s]');
ylabel('normalized Y_c [m]' );
saveas(gcf,[Main_Folder,'\','Y_movement.png']);


%---plot3 - marking of asperities by color
EventFrame = fixed_frames(:,:,EventFrame_num); 
EventFrame = [EventFrame;ones(size(EventFrame))];
figure; imshow(EventFrame);
hold on;
for asperityIdx = 1:length(meanVecX)
    plot(meanVecX(asperityIdx),12,'o',...
        'MarkerFaceColor',ColorSet(asperityIdx,:),...
        'Color',ColorSet(asperityIdx,:),...
        'MarkerSize',7);
end

% %---plot4 - superimposing the asperities trajectories for easy comparison
% figure;
% hold on;
% for i = 1:length(CoordinateCell)
%     mean_pixel_X = mean(Xcenter(i,1:pre_frames));
%     mean_pixel_Y=mean(Ycenter(i,1:pre_frames));
%     
%     plot3(Xcenter(i,:)-mean_pixel_X,...
%         Ycenter(i,:)-mean_pixel_Y,...
%         i*10*ones(size(Xcenter(i,:))),...
%         'Color',ColorSet(i,:),...
%         'LineWidth',3);
% end
% 
% %---plot 5 - graphing speeds
% XX=diff(Xcenter,1,2);           % find speeds
% MaxSpeedXX = max(abs(XX),[],2); % maximum speed magnitude
% figure;
% hold all;
% for i=1:length(MaxSpeedXX)
% plot(i,MaxSpeedXX(i),...
%         'o','Color',ColorSet(i,:),'MarkerFaceColor',ColorSet(i,:));
% end
% plot(1:length(MaxSpeedXX),MaxSpeedXX);

%----plot6 - mean movement on X direction
 
figure;
plot(Time,mean(Xcenter_normalized));
title('mean movement on X direction');
xlabel('Time [s]');
ylabel('Position [m]');
ylim([-3e-6 3e-6]);
% Vmax=max(diff(mean(Xcenter_normalized)));

figure;
plot(Time,mean(Ycenter_normalized));
title('mean movement on Y direction');
xlabel('Time [s]');
ylabel('Position [m]');
ylim([-3e-6 3e-6]);
% Vmax=max(diff(mean(Xcenter_normalized)));

%% save additional data
%--- add data to movie data structure
% analyzeDataStruct.MaxSpeedsX = MaxSpeedXX;
analyzeDataStruct.meanVecX = meanVecX;
analyzeDataStruct.meanVecY = meanVecY;

%--- overwrite the data structure on the computer:
save(video_struct, 'analyzeDataStruct');

%%  mass_center plots
% %plot1 -  Xcenter
% figure
% hold all
% for regionIdx=1:length(coordinatesCell)
%     plot(frame_range,Xcenter(regionIdx,:),'.')
% end
% hold off
% title('\fontsize{16}Mass-Center location on X axis');
% xlabel('Frames');
% ylabel('X_c [pixel]' );

% %%plot2 -  Ycenter
% figure
% hold all
% for regionIdx=1:length(coordinatesCell)
%     plot(frame_range,Ycenter(regionIdx,:),'.')
% end
% hold off
% title('\fontsize{16}Mass-Center location on Y axis');
% xlabel('Frames');
% ylabel('Y_c [pixel]' );

% %%plot3 - scatter
% figure; 
% hold on;
% for regionIdx = 1:length(coordinatesCell)
%     scatter(Xcenter(regionIdx,1:end-1),Ycenter(regionIdx,1:end-1),[],frame_range(1:end-1),'filled');          %color is time
% end
% caxis([frame_range(1),frame_range(end)]);       %colorbar for time
% colorbar;

%% peak plots 

%%% plot1 - scatter peaks
% TraceFigure = figure; 
% hold on;
% for regionIdx = 1:length(coordinatesCell)
% %     scatter(MaxXcoor(regionIdx,1:end-1), MaxYcoor(regionIdx,1:end-1),[],speedMag(regionIdx,:),'filled');        %color is speed
%     scatter(MaxXcoor(regionIdx,1:end-1), MaxYcoor(regionIdx,1:end-1),[],frame_range(1:end-1),'filled');          %color is time
% end
% % caxis([min(speedMag(:)),max(speedMag(:))]);   %colorbar for speed
% caxis([frame_range(1),frame_range(end)]);       %colorbar for time
% colorbar;


%%% plot2 - Peak location on Y axis
% figure
% hold all
% for i=1:6
% plot(frame_range,MaxYcoor(i,:))
% end
% title('\fontsize{16}Peak location on Y axis');
% xlabel('Frames');
% ylabel('Ypixel' );

%%% plot3 - Peak location on X axis
% figure
% hold all
% for regionIdx=1:length(coordinatesCell)
% plot(frame_range,MaxXcoor(regionIdx,:),'.')
% end
% title('\fontsize{16}Peak location on X axis');
% xlabel('Frames');
% ylabel('Xpixel' );

%% speed plots

% %%plot - total speed
% figure
% hold all
% for i=1:6
% plot(frame_range(1:80),speedMag(i,:))
% end
% title('\fontsize{16}Total Speed');
% xlabel('Frames');
% ylabel('Ypixel' );
