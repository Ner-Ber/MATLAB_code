function [RowOverTimeCell, vallyPairsCell,asper_num, RowOverTimeNormCell,...
    maxTrajectories, CM_Trajectories, skewnessCell, varCell,...
    center_of_massCell, sg_data, Images,Cam_Meta]...
    = analyzeEnlargementAndStrain(experNum, eventNum,varargin)
%analyzeEnlargementAndStrain will load and plot data from camera and strain
% and follow asperities from movies

%% define defaults
[pre_time, post_time, lineNum] = setDefaults4function(varargin,0.5*1e-3, 4e-3,'all');

exper=my_dir;
disp('');

%% Create the Images matrix, [rows * spatial * time]
disp('Reading images...');
folderInExperDir = my_dir(exper{experNum});
Cam_Meta = CameraMetaAllCams(exper{experNum});
expDetails=expDetailsRead(exper{experNum});
triggerDelay=expDetails.triggerDelay;
pre_frames = round(pre_time*Cam_Meta.FrameRate);
post_frames = round(post_time*Cam_Meta.FrameRate);
if ismember('Ph',folderInExperDir)
    eventFrame = Cam_Meta.NumIms-(Cam_Meta.PostIms+round(triggerDelay*Cam_Meta.FrameRate));
    if (eventFrame-pre_frames)<1;    fromFrame=1; else	fromFrame=eventFrame-pre_frames;	end
    if post_frames>Cam_Meta.PostIms;	toFrame=Cam_Meta.NumIms;	else	toFrame=eventFrame+post_frames;	end
    Images = phantomReadIms(exper{experNum},eventNum,fromFrame,1,toFrame,1,lineNum);
%     Images = phantomReadIms(exper{experNum},eventNum,1,1,Cam_Meta.NumIms,1,lineNum);
%     Images = Images(:,:,fromFrame:toFrame);
elseif ismember('IDT',folderInExperDir)
    [Images,~] = IDT_ReadImages(exper{experNum}, eventNum,pre_frames,1,post_frames,0,lineNum);   %IDT
else
    [Images,~] = photronReadImages(exper{experNum},eventNum,pre_frames,1,post_frames,0,lineNum);   %Photron
    
end
%--- average on rows:
Images = mean(Images,1);
flip_images =0;
if flip_images ==1
    Images = -Images;
    Images = Images+min(Images(:));
end

%% Make The Row-Over-Time Pictures, each cell for a different row
disp('creating RowOverTimeCell...');
Rows=1:size(Images,1);                 %% choose rows to look at
RowOverTimeCell = cell(1,length(Rows));
RowOverTimeCell{1} = permute(Images,[3 2 1]);

% Mean_RowOverTime=mean(Images,1);
% Mean_RowOverTime=permute(Mean_RowOverTime,[3 2 1]);
% for rowIdx = 1:length(Rows)
%     RowOverTime = Images(rowIdx,:,:);
%     RowOverTime = permute(RowOverTime,[3 2 1]);
%     RowOverTimeCell{rowIdx} = RowOverTime;
% end
% RowOverTimeCell{length(Rows)+1}=Mean_RowOverTime;
% Rows = [Rows, Rows(end)+1];
% Images = cat(1,Images,mean(Images,1));

%% Define the Asperities to follow

disp('Definning and following asperities...');
% factor=1.3;            % factor for threshold, asperities will be choosed only if they are higher than it
% min_gap=5;             % minimum gap between asperities
% asperity_min_size=12;  % minimum width of asperities
% [vallyPairsCell,asper_num]=Movie_define_asperities(RowOverTimeCell,Rows,factor,min_gap,asperity_min_size);

[vallyPairsCell,asper_num] = Movie_define_asperities_for_scratches(RowOverTimeCell);
disp('num of asperities in each row:');
disp(asper_num);

% phedislocationInPix = Movie_findTrenches(relevanrRowOverTime(270,:));

% vallyPairsCell = my_vallyPairsCell;

%% Make the images ready for analizing

MinOrMax = 'min';       % max to follow trench, 'min' to follow ridge
[RowOverTimeNormCell,~, ~,MarkedMatCell, ~] =...
    Movie_follow_maxima_in_row(Images,vallyPairsCell, Rows, 'all', [], 1, [],MinOrMax);   % all other cells

%% Analyze the picture to get the trajectories of desired data:
[maxTrajectories, CM_Trajectories, skewnessCell, varCell,center_of_massCell] =...
    Movie_analyze_asperities(MarkedMatCell, RowOverTimeCell);

%% Get data from SG:
disp('Reading strain gages data...');
sg_data=acq132_event_get_data...
    (exper{experNum},eventNum,'start','end',1,'Uxx','Uyy','Uxy','x_sg','y_sg','F','N');

end