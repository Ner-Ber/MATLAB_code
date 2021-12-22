function [minLag,RMS,shiftVector] = xcorr_vertical(x,y,shiftParamter)
%[r,lags] = XCORR_VERTICAL(x,y,varargin)
%
% XCORR_VERTICAL will shift the signal y in the vertical direction to find
% the shift which gives the maximal correlation.
% shiftParamter - shift steps can be set be entering a shifts-vector or by
% entering a scalar, N, which means that the shifts examined will be of
%  shiftVector=[(-min(y)+max(x) : (peak2peak(x)/N) : (-max(y)+min(x))]


%% create shifts vector
if nargin<3
    ds = (peak2peak(x)/1000);
    shiftVector=((-max(y)+min(x)-ds):ds:(-min(y)+max(x))+ds);
else
    if isscalar(shiftParamter)
        ds = (peak2peak(x)/shiftParamter);
        shiftVector=((-max(y)+min(x)-ds):ds:(-min(y)+max(x))+ds);
    elseif isvector(shiftParamter)
        shiftVector = shiftParamter;
    else
        error('invalid third input')
    end
end

%% find correlation between shifts
%--- create matric containing all shifts of y:
y_shifts = bsxfun(@plus,repmat(y(:)',length(shiftVector),1),shiftVector(:));
DiffXY = bsxfun(@minus,y_shifts,x(:)');
%--- omit nans:
DiffXY(isnan(DiffXY)) = 0;
%--- calc RMS:
RMS = rms(DiffXY,2);
%--- find minimizig lag
[~,I] = min(RMS);
minLag = shiftVector(I);
end


