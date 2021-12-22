function imsStruct = IDT_imagesWithData(exp_dir,eventNum,pre_time,intervalTime,post_time,smt,lineNum)
%imsStruct = IDT_imagesWithData(exp_dir,eventNum,pre_time,intervalTime,post_time,smt,lineNum)
% will create a structure with lines (row over time) and more data.
%exp_dir - path of the experiment folder
%eventNum - event number
%pre_time - time in seconds to shoe before event
%intervalTime - time in seconds between frames
%post_time - time in seconds to show after event
%smt - smooth factor along the x axis (default=0)
%lineNum - lines to consider (default = 'all');

%--------Read Phantom Meta
meat_path= [exp_dir '\IDT'];
IDT_Meta = IDT_ReadMeta(meat_path);

%--- set frames from pre and post times
pre_event = round(pre_time/(IDT_Meta.FrameT*1e-6));
if pre_event>(IDT_Meta.NumIms-IDT_Meta.PostIms)
    pre_event=(IDT_Meta.NumIms-IDT_Meta.PostIms);
end

interval = round(intervalTime/(IDT_Meta.FrameT*1e-6));
if (interval>IDT_Meta.NumIms) || interval<1
    interval=1;
end

post_event = round(post_time/(IDT_Meta.FrameT*1e-6));
if post_event>(IDT_Meta.PostIms)
    post_event=IDT_Meta.PostIms;
end

%--- read images
[Images,CAM_meta] = IDT_ReadImages(exp_dir,eventNum,pre_event,interval,post_event,smt,lineNum);

%% build the structure
path=pwd;
lastSlashPosition = find(path == '\', 1, 'last');
imsStruct.Date = path(lastSlashPosition+1:end);

imsStruct.Exp=exp_dir;
imsStruct.EventNum=eventNum;

imsStruct.trigger='NaN';
imsStruct.t = (-pre_event:interval:post_event)*IDT_Meta.FrameT*1e-6;    % in seconds
imsStruct.lines = permute(mean(Images,1),[3 2 1]);
imsStruct.included = lineNum;

exp_details=expDetailsRead(exp_dir);
expMeasures = expMeasureProperties(exp_dir);
res = expMeasures.resBigPic;
blockEdgesCorrdinates = expMeasures.blockEdgesCorrdinates;
% blockEdgesCorrdinates = IDT_findEdgePixlesOfBlock(imsStruct.lines);
% blockEdgesCorrdinates = [0 1920];
% xx = linspace(0,exp_details.UpperBlockLength*1e-3,abs(diff(blockEdgesCorrdinates)));
imsStruct.res = res; % meters per pixel
alongFrame = 1:IDT_Meta.ImageWidth;
alongFrameFixed = alongFrame-blockEdgesCorrdinates(1);
alongFrameRes = alongFrameFixed*imsStruct.res;
imsStruct.x = alongFrameRes; % in meters
% imsStruct.x = linspace(0,exp_details.UpperBlockLength*1e-3,IDT_Meta.ImageWidth); % in meters
imsStruct.fps = CAM_meta.FrameRate;
imsStruct.BlockEdgePixles = blockEdgesCorrdinates;

end