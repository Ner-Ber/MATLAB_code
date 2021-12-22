function [T_mean,PhediLoc_mean] = phedi_averagePhedi(DataStruct, phediIdxs, varargin)
%[avg_phedi] = phedi_averagePhedi(dataStruct, phediIdxs, reduceSin, t_tips, plotSinFit)

%% set defaults
[reduceSin, t_tips, plotSinFit] = setDefaults4function(varargin,1,1,0);

%% create fixed matrix
PhediLoc = DataStruct.PhediData.PhediLocation;
if t_tips
    T = bsxfun(@minus,DataStruct.PhediData.timeVec,DataStruct.PhediData.t_tips(:)');
else
    T = DataStruct.PhediData.timeVec;
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

%% fit shifted figures;

%--- interpolate for common time steps
t_unique = unique(T);
T_phedis = repmat(t_unique(:),1,N);
PhediLocFixed_iterped = nan(size(T_phedis));
for i=1:N
    PhediLocFixed_iterped(:,i) = interp1(T(:,i),PhediLocFixed(:,i),T_phedis(:,i));
end
%--- create mean measurement
T_mean = mean(T_phedis(:,phediIdxs),2);
PhediLoc_mean = mean(PhediLocFixed_iterped(:,phediIdxs),2);

end