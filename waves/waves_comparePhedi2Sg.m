function waves_comparePhedi2Sg(DataStruct,varargin)
% waves_comparePhedi2Sg(DataStruct,strain_component,minTime,maxTime)
% This works only for non-events recordings

%% set defaults 
[strainComp,minTime,maxTime] = setDefaults4function(varargin,'Uxx',-inf,inf);

SgData= DataStruct.SgData;
t = SgData.t*1e-3;

if strcmpi(strainComp,'Uxx')
    StrainMat = SgData.Uxx;
elseif strcmpi(strainComp,'Uxy')
    StrainMat = SgData.Uxy;
elseif strcmpi(strainComp,'Uyy')
    StrainMat = SgData.Uyy;
end
%% process and prepare
timeVec = DataStruct.PhediData.timeVec;
%--- crop
minSg = min(t);
maxSg = max(t);
minPhedi = min(timeVec);
maxPhedi = max(timeVec);
begTime = max([minSg,minPhedi,minTime]);
endTime = min([maxSg,maxPhedi,maxTime]);
%--- sg corrdinated
[~,begSgIdx] = min(abs(t-begTime));
[~,endSgIdx] = min(abs(t-endTime));
[~,begPhIdx] = min(abs(timeVec-begTime));
[~,endPhIdx] = min(abs(timeVec-endTime));

tCrop = t(begSgIdx:endSgIdx);
StrainCrop = StrainMat(begSgIdx:endSgIdx,:);
StrainCrop = bsxfun(@minus,StrainCrop,mean(StrainCrop(1:70,:)));
%--- super smooth
strainCropSmthMat = [];
for i=1:15
    strainCropSmthMat(:,i) = supsmu(tCrop,StrainCrop(:,i));
end
%--- sum to get displacement 
SummedMat = cumsum(strainCropSmthMat,1)*1e-6;
%--- normalize
ShiftedMat = bsxfun(@minus,SummedMat,min(SummedMat,[],1));
StrchdMat = bsxfun(@rdivide,ShiftedMat,max(ShiftedMat,[],1));
StrchdMat = StrchdMat*2-1;

%% plot 
N_sg = size(StrchdMat,2);
SG_colors = MyVaryColor(N_sg,jet);
wavesDisplay_detrend_and_strech(DataStruct,[begPhIdx:endPhIdx],'RawNorm');
hold on;
for i=1:N_sg
    plot(tCrop,StrchdMat(:,i),'-','Color',SG_colors(i,:),'LineWidth',3);
end
