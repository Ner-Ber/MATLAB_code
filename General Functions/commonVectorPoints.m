function [commonX, Y1c, Y2c] = commonVectorPoints(X1,Y1,X2,Y2,varargin)
% [commonX, Y1c, Y2c] = commonVectorPoints(X1,Y1,X2,Y2)
%
% commonVectorPoints will create a common X axis of the unique point of
% data sets 1,2, and interpolate (linear) the data in Y1,2 to fit the
% common axis.

[res] = setDefaults4function(varargin,'none');

commonXunq = unique([X1(:);X2(:)]);
if ~ischar(res)
    commonXunq = unique(round(commonXunq,-floor(log10(res))));
end
maximalCoor = min(max(X1(:)),max(X2(:)));
minimalCoor = max(min(X1(:)),min(X2(:)));
commonX = commonXunq;
commonX(~(commonXunq>minimalCoor & commonXunq<maximalCoor)) = [];
Y1c = interp1(X1,Y1,commonX);
Y2c = interp1(X2,Y2,commonX);

end
