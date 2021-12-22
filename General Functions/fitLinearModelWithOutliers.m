function [slope, intersection, stdVec, iterations] = fitLinearModelWithOutliers(X,Y)
% [slope, intersection, stdVec, iterations] = fitLinearModelWithOutliers(X,Y)

newX = X;
newY = Y;
stdVec = [];

maxIter = 1e3;
iterations = 0;
converged = 0;

while ~converged && iterations<maxIter
    %--- fit model
    P = polyfit(newX,newY,1);
    modelY = P(2)+P(1)*newX;
    %--- get residuals
    residuals = abs(modelY-newY);
    STD = std(residuals);
    %--- find outliers
    outilerCutoff = 3*STD;
    NotOutliersLogical = residuals<outilerCutoff;
    %--- update process
    stdVec = cat(1,stdVec,STD);
    %--- update data vectors
    newX = newX(NotOutliersLogical);
    newY = newY(NotOutliersLogical);
    
    %--- check convergence
    STD_thresh = 1;
    converged = STD<STD_thresh;
    %--- propagate iteration counter 
    iterations = iterations+1;
    
end

intersection = P(2);
slope = P(1);


end