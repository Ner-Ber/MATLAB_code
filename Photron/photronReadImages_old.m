function [Images] = photronReadImages(exp_dir,eventNum,startl,interval,endl,smt,includeFlag,lineNum)
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
for i = startl:interval:endl
    %--- read image
    f = fopen([eventPath '\' ImagesList{i}],'r');
    Im=fread(f,imSize,'uint16');
    fclose(f);
    
    %--- reshape and edit if needed:
    Imr=reshape(Im',PhMeta.ImageWidth,PhMeta.ImageHeight)';
    Imr=Imr(lineNum,:);
    Imr = Imr/2^16;
    if isscalar(smt) && (smt>0)
        Imr = photronSmoothImages(Imr, smt);
    end
    Images = cat(3,Images,Imr);
end


