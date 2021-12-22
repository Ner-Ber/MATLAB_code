function acqOut=smooth_sg_x(acq,smtX)
%The function get acq132 struct and smoothes all data computed from sg in x
%direction

if nargin<2
    smtX=3;
end

if length(acq.t)>200
    indexEnd=200;
else
    indexEnd=length(acq.t);
end

fNames=fieldnames(acq);
for j=1: length(fNames)
    if (length(acq.(fNames{j})(:,1))==length(acq.t(:,1))) && (length(acq.(fNames{j})(1,:))==length(acq.x_sg))
        tmp=acq.(fNames{j});
        pre=mean(tmp(1:indexEnd,:),1);
        preSmt=smooth(pre,smtX)';
        acqOut.(fNames{j})=acq.(fNames{j})-repmat(pre,length(acq.t),1)+repmat(preSmt,length(acq.t),1);
    else
        acqOut.(fNames{j})=acq.(fNames{j});
    end
end



