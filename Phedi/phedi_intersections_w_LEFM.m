% phedi_intersections_w_LEFM(PhediStruct)
%
% phedi_intersections_w_LEFM will find the intersection of the mean phedi
% (displacement and particle velocity) with the LEFM prediction (which is
% based on Cf and fracture energy).
% Intended to find intersections in the negativ x<0, thus behind crack tip.
%
% VARARGIN = DAFAULTS:
% plotSlopes='11'       (two digits:
%                       1st-indicates slope to plot (+/-1) 0 means all slopes.
%                       2nd-logical indicating if plot damaged.

function [maxXi_loc,maxYi_loc,maxXi_vel,maxYi_vel] = phedi_intersections_w_LEFM(DataStruct,varargin)
    
    %% set parameters
    [plotSlopes] = setDefaults4function(varargin,'11');
    
    
    PhediData = DataStruct.PhediData;
    PhediSlopes = PhediData.slopeIncline;
    PhediSlopes = mean(PhediSlopes,1);
    PhediSlopes(~~mod(PhediSlopes,1)) = 0;
    plotDamaged = str2double(plotSlopes(end));
    slope2Plot = str2double(plotSlopes(1:end-1));
    phdIndx = [];
    if plotDamaged
        phdIndx = cat(2,phdIndx,find(PhediSlopes==0));
    end
    if slope2Plot==-1
        phdIndx = cat(2,phdIndx,find(PhediSlopes==-1));
    elseif slope2Plot==1
        phdIndx = cat(2,phdIndx,find(PhediSlopes==1));
    else
        phdIndx = cat(2,phdIndx,find(PhediSlopes~=0));
    end
    phdIndx = unique(phdIndx);
    
    
    
    AvgPhediStruct = phedi_averagePhedi(DataStruct,phdIndx);
    %% find intersections with particle velocity:
    negativeLogic_Avg = AvgPhediStruct.X_mean<0 & AvgPhediStruct.X_mean>-0.05;
    negativeLogic_LEFM = DataStruct.solAtInter.x<0 & DataStruct.solAtInter.x>-0.05;
    
    [Xi_vel,Yi_vel] = intersections(AvgPhediStruct.X_mean(negativeLogic_Avg),AvgPhediStruct.Vel_mean_x(negativeLogic_Avg),...
        DataStruct.solAtInter.x(negativeLogic_LEFM),DataStruct.solAtInter.vx(negativeLogic_LEFM));
    
    [maxXi_vel,I] = max(Xi_vel);
    maxYi_vel = Yi_vel(I);
    
    
    %% find intersections with particle location:
    
    [Xi_loc,Yi_loc] = intersections(AvgPhediStruct.X_mean(negativeLogic_Avg),AvgPhediStruct.Loc_mean_x(negativeLogic_Avg),...
        DataStruct.solAtInter.x(negativeLogic_LEFM),DataStruct.solAtInter.Ux(negativeLogic_LEFM));
    
    [maxXi_loc,I] = max(Xi_loc);
    maxYi_loc= Yi_loc(I);
    
    if isempty(maxXi_vel); maxXi_vel=nan; end;
    if isempty(maxYi_vel); maxYi_vel=nan; end;
    if isempty(maxXi_loc); maxXi_loc=nan; end;
    if isempty(maxYi_loc); maxYi_loc=nan; end;
end