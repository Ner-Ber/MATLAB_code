function [Mov4D_colorMat] = IDT_createNormalizedMovie(exp_dir,eventNum,varargin)

[pre_time,intervalTime,post_time,smt,lineNum] = setDefaults4function(varargin,...
        3e-3,1e-6,3e-3,0,'all');

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

%% read images
[Images,~] = IDT_ReadImages(exp_dir,eventNum,pre_event,interval,post_event,smt,lineNum);
%% reshape into 2D for colormapping
ImagesNormed = bsxfun(@rdivide,Images,Images(:,:,1));
S = size(Images);

ImagesLin = reshape(ImagesNormed,S(1)*S(2),[]);

%% normalize for colormapping
cMax = 1.2;
cMin = 0;
ImagesLin(ImagesLin<cMin)=cMin;
ImagesLin(ImagesLin>cMax)=cMax;
ImagesLin = ImagesLin./cMax;
ImagesLin = uint8(floor(ImagesLin*255));
%% do colormapping
Map = MyVaryColor(255,jet);
ColoeredLin = ind2rgb(ImagesLin, Map);
ColoeredLin = permute(ColoeredLin,[1 3 2]);
%% reshape
Mov4D_colorMat = reshape(ColoeredLin,S(1),S(2),3,S(3));

end