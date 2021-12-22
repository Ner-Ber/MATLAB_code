function DataStructFixed = phedi_reduceSineFromPhediTraj(DataStruct,varargin)

[reduceSin, t_tips, plotSinFit] = setDefaults4function(varargin,1,1,0);


PhediLoc = DataStruct.PhediData.PhediLocation;
if t_tips
    T = bsxfun(@minus,DataStruct.PhediData.timeVec,DataStruct.PhediData.t_tips(:)');
else
    T = PhediData.timeVec;
end

N = size(T,2);
if reduceSin
    myfitStruct = phedi_fitSineToPhediTraj(DataStruct.PhediData,[],2,plotSinFit);
    PhediLocFixed = nan(size(PhediLoc));
    for i=1:N
        PhediLocFixed(:,i) = PhediLoc(:,i)-myfitStruct.model(myfitStruct.b_param{i},myfitStruct.Shifts(i),T(:,i));
    end
else
    PhediLocFixed = PhediLoc;
end


DataStructFixed = DataStruct;
DataStructFixed.PhediDataOriginal = DataStruct.PhediData;
DataStructFixed.PhediData.PhediLocation = PhediLocFixed;