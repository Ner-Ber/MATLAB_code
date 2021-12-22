function updated_sgDataStruct = sg_createVectorsForSg(sgDataStruct,BigPicRotStruct)
%sg_createVectorsForSg will update the structure containing the strain-gage
%data so it will contain info about the from passing at each sg.
%
% INPUTS:
% sgDataStruct - is the output of the function acq132_event_get_data
% BigPicRotStruct - is the output of the fucntion IDT_PlotRowOverTime
%
%The function will add the fields:
%v_sg = velocity of front at sg location
%t_tips_sg = time front passed the sg location;
%t_mins_t_tips = time vecotr for sg of t-t_tip;
%x_mins_x_tips = location vecotr for sg of x-x_tip;

%% find velocity at locations of sg
%--- create a location vector for velocity by mapping from location to pixles
LocationVelocityVector = interp1(1:length(BigPicRotStruct.x),BigPicRotStruct.x,BigPicRotStruct.frontVelLoc);
velocityAtSg = interp1(LocationVelocityVector,BigPicRotStruct.frontVel*BigPicRotStruct.fps*BigPicRotStruct.res,sgDataStruct.x_sg/1000);

%% create timevectors of t-t_tip (this assumes time in phantom and in peter are sync)
N = length(sgDataStruct.x_sg);
timeMat = repmat(sgDataStruct.t(:)/1000,1,N);
[~,sg_indexes_4Loc] = min(abs(bsxfun(@minus,repmat(BigPicRotStruct.x(:),1,N),sgDataStruct.x_sg/1000)));
t_tips = BigPicRotStruct.frontTime_interp(sg_indexes_4Loc)/BigPicRotStruct.fps;
t_mins_t_tips = bsxfun(@minus, timeMat,t_tips(:)');

%% create x-x_tip vectors
x_mins_x_tips = -bsxfun(@times,t_mins_t_tips,velocityAtSg);

%% update structure
updated_sgDataStruct = sgDataStruct;
updated_sgDataStruct.v_sg = velocityAtSg;
updated_sgDataStruct.t_tips_sg = t_tips;
updated_sgDataStruct.t_mins_t_tips = t_mins_t_tips;
updated_sgDataStruct.x_mins_x_tips = x_mins_x_tips;