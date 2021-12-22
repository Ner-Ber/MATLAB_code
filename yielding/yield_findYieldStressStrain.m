%S = yield_findYieldStressStrain(DataStruct,varargin)
%[averagingLength,phediIdxs] = setDefaults4function(varargin,1e-3,true(1,N));
%

function S = yield_findYieldStressStrain(DataStruct,varargin)
    %% set defaults
    PhediVel = DataStruct.PhediData.PhediVelocity;
    N = size(PhediVel,2);
    
    [averagingLength,phediIdxs] = setDefaults4function(varargin,1e-3,true(1,N));
    %% get properties of velocity on interface
    AvgPhediStruct  = phedi_averagePhedi(DataStruct,phediIdxs);
    S = phedi_findRisingPointFromVel(AvgPhediStruct.Vel_mean_x,AvgPhediStruct.X_mean);
    S.frstSlope = S.VelKinkY./(S.VelKinkX-S.risingVelPoint);
    S.scndSlope= (S.v0-S.VelKinkY)/(-S.VelKinkX);
    
    %% set smooth
    smthL = averagingLength/abs(diff(DataStruct.solAtInter.x([1,2])));
    
    %% get properties of stress and strain from LEFM on interface
    [~,S.yieldSxx1st] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Sxx,smthL),S.risingVelPoint*[1 1],inf*[-1 1]);
    [~,S.yieldSyy1st] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Syy,smthL),S.risingVelPoint*[1 1],inf*[-1 1]);
    [~,S.yieldSxy1st] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Sxy,smthL),S.risingVelPoint*[1 1],inf*[-1 1]);
    [~,S.yieldUxx1st] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Uxx,smthL),S.risingVelPoint*[1 1],inf*[-1 1]);
    [~,S.yieldUyy1st] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Uyy,smthL),S.risingVelPoint*[1 1],inf*[-1 1]);
    [~,S.yieldUxy1st] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Uxy,smthL),S.risingVelPoint*[1 1],inf*[-1 1]);
    
    [~,S.yieldSxx2nd] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Sxx,smthL),S.VelKinkX*[1 1],inf*[-1 1]);
    [~,S.yieldSyy2nd] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Syy,smthL),S.VelKinkX*[1 1],inf*[-1 1]);
    [~,S.yieldSxy2nd] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Sxy,smthL),S.VelKinkX*[1 1],inf*[-1 1]);
    [~,S.yieldUxx2nd] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Uxx,smthL),S.VelKinkX*[1 1],inf*[-1 1]);
    [~,S.yieldUyy2nd] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Uyy,smthL),S.VelKinkX*[1 1],inf*[-1 1]);
    [~,S.yieldUxy2nd] = intersections(DataStruct.solAtInter.x,smooth(DataStruct.solAtInter.Uxy,smthL),S.VelKinkX*[1 1],inf*[-1 1]);
    
    S.Cf = DataStruct.PhediData.Cf;
    
end
