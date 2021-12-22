function PhediLocReduced = phedi_reduceSinFromLocation(PhediData,varargin)
% PhediLocReduced = phedi_reduceSinFromLocation(PhediData, numOfSin, preKinkTime, DoPlot)


[numOfSin, preKinkTime, DoPlot] = setDefaults4function(varargin,2,-0.025/PhediData.Cf,0);

PhediLoc = PhediData.PhediLocation;
T = bsxfun(@minus,PhediData.timeVec(:),PhediData.t_tips(:)');

N = size(T,2);
myfitStruct = phedi_fitSineToPhediTraj(PhediData,preKinkTime,numOfSin,DoPlot);
PhediLocReduced = nan(size(PhediLoc));
for i=1:N
    PhediLocReduced(:,i) = PhediLoc(:,i)-myfitStruct.model(myfitStruct.b_param{i},myfitStruct.Shifts(i),T(:,i));
end


end