function Rsq = r_square(yData,yModel)
% Rsq = r_square(yData,yModel)

SStot = sumsqr(yData-mean(yData));
SSres = sumsqr(yData-yModel);
Rsq = 1-(SSres/SStot);

end