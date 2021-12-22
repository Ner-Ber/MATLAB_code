function [Images,IDT_Meta] = IDT_ReadImages(exp_dir,eventNum,pre_event,interval,post_event,smt,includeFlag,lineNum)
% [Images] = photronReadImages(exp_dir,eventNum,startl,interval,endl,smt,includeFlag,lineNum)
%
%eventNum should be a scalar.
%lineNum -> the vector of numbers of rows from the image.
%includeFlag -> 'all' for include all.
%returns [Lines Included]
%smt - avergae length for smoothing image. can be a two element vector
%which implies on the x and y legnth of the convoluting filter. 0 will
%indicate no smoothing. 

path= [exp_dir '\IDT'];
% path=[path '\TestSession_010'];   %***just for the 23/8 experiment***

%--------Read Phantom Meta
IDT_Meta = IDT_ReadMeta(path);
% IDT_Meta.BROCLenght=0;            %***just for the 23/8 experiment***

if (nargin<8 || strcmp(lineNum,'all')) %the full image is taken
    lineNum=(1:IDT_Meta.ImageHeight);
end

if(nargin<7 )
    includeFlag='';
end


%--- create list of all photos
ImagesList = struct2cell(dir([path,'/*.tif']));
ImagesList = ImagesList(1,:)';



%% determine event frame, starting frame and ending frame
Start_frame=IDT_Meta.StartFrame;
Frame_rate=IDT_Meta.FrameRate;
trig_delay=7e-3;                       % trig delay is 7 msec
event_frame=round(Start_frame-trig_delay*Frame_rate);
IDT_Meta.EventFrame=event_frame;

startl=event_frame-pre_event;
endl=event_frame+post_event;


%% ----- read images in to Images aray:
Images = [];
for i = ((eventNum-1)*IDT_Meta.BROCLenght+startl):interval:((eventNum-1)*IDT_Meta.BROCLenght+endl)
    %--- read image
    Imr = imread([path,'/',ImagesList{i}]);
    Imr=Imr(lineNum,:);
    %--- normalize to 1:
    Imr = double(Imr)/2^(2^ceil(log2(IDT_Meta.EffectiveBitDepth)));
    if isscalar(smt) && (smt>0)
        Imr = photronSmoothImages(Imr, smt);
    end
    Images = cat(3,Images,Imr);
end



end