function RowOverTime = phantom_getRowOverTime_BigPic(expPath, eventNum, Cam_Meta,varargin)
%% define defaults
[pre_time, post_time, lineNum] = setDefaults4function(varargin,0.5*1e-3, 4e-3,'all');

%% Create the Images matrix, [rows * spatial * time]
disp('Reading images...');
expDetails=expDetailsRead(expPath);
triggerDelay=expDetails.triggerDelay;
pre_frames = round(pre_time*Cam_Meta.FrameRate);
post_frames = round(post_time*Cam_Meta.FrameRate);
%--- read images from phantom
eventFrame = Cam_Meta.NumIms-(Cam_Meta.PostIms+round(triggerDelay*Cam_Meta.FrameRate));
%--- trim extended time requests
if (eventFrame-pre_frames)<1;    fromFrame=1; else	fromFrame=eventFrame-pre_frames;	end
if post_frames>Cam_Meta.PostIms;	toFrame=Cam_Meta.NumIms;	else	toFrame=eventFrame+post_frames;	end
Images = phantomReadIms_BigPic(expPath,eventNum,fromFrame,1,toFrame,1,lineNum);

%--- average on rows:
Images = mean(Images,1);

%--- flip if necessary
flip_images =0;
if flip_images ==1
    Images = -Images;
    Images = Images+min(Images(:));
end

RowOverTime = permute(Images,[3 2 1]);

end