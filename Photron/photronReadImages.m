function [Images,Cam_Meta] = photronReadImages(exp_dir,eventNum,pre_event,interval,post_event,smt,lineNum)
% [Images] = photronReadImages(exp_dir,eventNum,startl,interval,endl,smt,includeFlag,lineNum)
%
%eventNum if a scalar - defines the folder containing frames when ordered
%by name (small to large). if a string should name a specific event folder.
%lineNum -> the vector of numbers of rows from the image.
%includeFlag -> 'all' for include all.
%returns [Lines Included]
%smt - avergae length for smoothing image. can be a two element vector
%which implies on the x and y legnth of the convoluting filter. 0 will
%indicate no smoothing. 

path= [exp_dir '\Photron'];
%--------Read Phantom Meta
PhMeta = photronReadMeta(path);

if (nargin<8 || strcmp(lineNum,'all')) %the full image is taken
    lineNum=(1:PhMeta.ImageHeight);
end

if(nargin<7 )
    includeFlag='';
end

%----- create name list of events
content_names = dir_names(path);
event_names = content_names;
%------specify specific event
if ~ischar(eventNum)
    eventName = event_names{eventNum};
else
    eventName = eventNum;
end
eventPath = [path '\' eventName];


imSize=PhMeta.ImageHeight*PhMeta.ImageWidth;
% Images=zeros(PhMeta.ImageHeight, PhMeta.ImageWidth, floor((endl-startl)/interval)+1);   % define array to cantain images

%% ----- read images in to Images aray:
ImagesList = dir_names(eventPath);
%----- remove .cih file from the images list
ImagesList = ImagesList(cellfun(@isempty, strfind(ImagesList,'cih')));
Images = [];

%% determine event frame, starting frame and ending frame
% Start_frame=0;
% Frame_rate=PhMeta.FrameRate;
% trig_delay=7e-3;                       % trig delay is 7 msec
% event_frame=round(Start_frame-trig_delay*Frame_rate)+abs(PhMeta.StartFrame);
% PhMeta.EventFrame=event_frame;

startl=PhMeta.EventFrame-pre_event;
endl=PhMeta.EventFrame+post_event;

for i = startl:interval:endl
    %------for raw files-----%
    %--- read image
%     f = fopen([eventPath '\' ImagesList{i}],'r');
%     Im=fread(f,imSize,'uint16');
%     fclose(f);
%     
%     %--- reshape and edit if needed:
%     Imr=reshape(Im',PhMeta.ImageWidth,PhMeta.ImageHeight)';
%     Imr=Imr(lineNum,:);
%     Imr = Imr/2^16;
    
    %------for tif files-----%
    Imr = imread([eventPath '\' ImagesList{i}]);
    Imr=Imr(lineNum,:);
    Imr = double(Imr)/2^(2^ceil(log2(PhMeta.EffectiveBitDepth)));
    
    
    
    if isscalar(smt) && (smt>0)
        Imr = photronSmoothImages(Imr, smt);
    end
    Images = cat(3,Images,Imr);
end

Cam_Meta=PhMeta;
 

