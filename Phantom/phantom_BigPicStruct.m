function RowOverTimeStruct = phantom_BigPicStruct(exp_dir, eventNum, varargin)
% RowOverTimeStruct = phantom_BigPicStruct(exp_dir, eventNum, pre_time,intervalTime,post_time,timeBase,smt,lineNum, FallingPercnt)
%
% phantom_BigPicStruct will create a structure containing RowOverTime from
% the big optical path when using the phantom. 
% 
% INPUTS [DEFAULT VAL.]:
% exp_dir - local path of the wanted experiment
% eventNum - ...
% pre_time - time prior to trigger to account for (in absolute value) [0.5e-3]
% intervalTime - time between frames. usually want it as small as possible [1e-8]
% post_time - time prior to trigger to account for (in absolute value) [4e-3]
% timeBase - 1 for second [1]
% smt - length of linear smoothing mask [0, no smooth]
% lineNum - which line of frame to use ['all']
%

%% define defaults
% [pre_time, post_time, lineNum] = setDefaults4function(varargin,0.5*1e-3, 4e-3,'all');
[pre_time,intervalTime,post_time,timeBase,smt,lineNum, FallingPercnt] = setDefaults4function(varargin,...
    0.5e-3,1e-8,4e-3,1,0,'all',2);
%% get general info
exp_details = expDetailsRead(exp_dir);
CamMetaStructBig = CameraMetaAllCams(exp_dir,'PhBig');
%% read and create data
%--- get ims struct
ims = phantomGetLines_BigPic(exp_dir,eventNum,-pre_time,intervalTime,post_time,timeBase,smt,'all',lineNum);
%--- create ROT:
% RowOverTime = mean(ims.lines,1);
RowOverTime = ims.lines;
%--- normalized ROT:
NormalizedRowOverTime = bsxfun(@rdivide,RowOverTime,mean(RowOverTime(1:100,:)));
%--- calculate resolution
expMeasures = expMeasureProperties(exp_dir);
res = expMeasures.resBigPic;
blockEdgesCorrdinates = expMeasures.blockEdgesCorrdinates;
% blockEdgesCorrdinates = IDT_findEdgePixlesAllEvents(exp_dir);
% blockEdgesCorrdinates = IDT_findEdgePixlesOfBlock(RowOverTime);
% blockEdgesCorrdinates = [0 size(RowOverTime,2)];
% xx = linspace(0,exp_details.UpperBlockLength*1e-3,abs(diff(blockEdgesCorrdinates)));
% res = abs(xx(2)-xx(1)); % meters per pixel
%--- spacial vector
alongFrameFixed = (1:size(RowOverTime,2))-blockEdgesCorrdinates(1);
x = alongFrameFixed*res;
%--- find front and interpolate
% frontStructure = IDT_FindCfFromRowOverTime(ims.lines);
frontStructure = IDT_FindCfFromRowOverTimeControlFallingPercnt(ims.lines,[],FallingPercnt);
frontIntepStruct = IDT_interpolateFrontBetweenSteps(frontStructure);
%--- get velocity related vecs in [m/s] and [m]
frontVelLoc_interpM = interp1(1:length(frontStructure.fronLoc),x,frontIntepStruct.VelLoc,'linear','extrap');
frontVel_interpMperS = frontIntepStruct.FrontVel*CamMetaStructBig.FrameRate*res;
%--- find average velocity in scratch region
scratchRegionMeters = [0.07 0.09];
xq = linspace(scratchRegionMeters(1),scratchRegionMeters(2),50);
scartchLocationsInPix = interp1(x,1:length(x),xq);
scartchLocationsInPix = unique([scartchLocationsInPix,round(min(scartchLocationsInPix)):round(max(scartchLocationsInPix))]);
frontTimeInterped = interp1(frontStructure.StepsPix,frontStructure.StepTime,scartchLocationsInPix);
AvgCf_scratches = mean(diff(x))*mean(diff(scartchLocationsInPix)./diff(frontTimeInterped))./mean(diff(ims.t));



% RowOverTime = phantom_getRowOverTime_BigPic(exp_dir, eventNum, CamMetaStructBig, Param.preRowsTime, Param.postRowsTime, Param.lineNum);


%% save to struct
RowOverTimeStruct = struct;
RowOverTimeStruct.DataMat = RowOverTime;
RowOverTimeStruct.DataMatNorm = NormalizedRowOverTime;
RowOverTimeStruct.res = res;
RowOverTimeStruct.fps = CamMetaStructBig.FrameRate;
RowOverTimeStruct.x = x;
RowOverTimeStruct.t = ims.t;
RowOverTimeStruct.details = [ims.Date,'   ',ims.Exp,'   event=',num2str(ims.EventNum)];
RowOverTimeStruct.xlabel = 'x [m]';
RowOverTimeStruct.ylabel = 't [s]';
RowOverTimeStruct.frontLoc = frontStructure.fronLoc;
RowOverTimeStruct.frontStepsPix = frontStructure.StepsPix;
RowOverTimeStruct.frontStepsTime = frontStructure.StepTime;
RowOverTimeStruct.frontVel = frontStructure.FrontVel;
RowOverTimeStruct.frontVelLoc = frontStructure.VelLoc;
RowOverTimeStruct.frontLoc_interp = frontIntepStruct.xq;
%     RowOverTimeStruct.frontTime_interp = frontIntepStruct.tq;
RowOverTimeStruct.frontTime_interp = frontIntepStruct.tq+ims.t(1)*CamMetaStructBig.FrameRate-1;    % time vector for front in units of frames
RowOverTimeStruct.frontTimeOriginal_interp = frontIntepStruct.tq_original+ims.t(1)*CamMetaStructBig.FrameRate-1;    % same as frontTime_interp but not smoothed
RowOverTimeStruct.frontVel_interp = frontIntepStruct.FrontVel;
RowOverTimeStruct.frontVelLoc_interp = frontIntepStruct.VelLoc;
RowOverTimeStruct.frontVelLoc_interpM = frontVelLoc_interpM;
RowOverTimeStruct.frontVel_interpMperS = frontVel_interpMperS;
RowOverTimeStruct.AvgCf_scratches = AvgCf_scratches;
RowOverTimeStruct.CamMetaStructBig = CamMetaStructBig;
RowOverTimeStruct.BlockEdgePixles = blockEdgesCorrdinates;
RowOverTimeStruct.FallingPercnt = FallingPercnt;
end