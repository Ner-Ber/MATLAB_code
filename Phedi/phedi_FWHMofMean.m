function [FWHM,FWHM_Y,Loc2Right,Loc2Left] = phedi_FWHMofMean(DataStruct,varargin)
    
    %% defaults
    PhediVel = DataStruct.PhediData.PhediVelocity;
    N = size(PhediVel,2);
    
    [averagingLength, phediIdxs] = setDefaults4function(varargin,1e-3,true(1,N));
    
    
    %%
    AvgPhediStruct  = phedi_averagePhedi(DataStruct, phediIdxs);
    relevantLogic = abs(AvgPhediStruct.X_mean)<=0.04;
    Vel_mean_xRelevant = AvgPhediStruct.Vel_mean_x(relevantLogic);
    X_meanRelevant = AvgPhediStruct.X_mean(relevantLogic);
    
    
    Vel_mean_xSmth = smooth(Vel_mean_xRelevant,max(1,round(averagingLength/abs(diff(X_meanRelevant(1:2))))));
    
    %-- find peak:
    [~,I] = max(Vel_mean_xSmth);
    peakLoc = X_meanRelevant(I);
    peakval = Vel_mean_xRelevant(I);
    %--- fwhm val:
    fwhmVal = peakval/2;
    %--- fwhm loc:
    [X0, ~] = intersections(X_meanRelevant,Vel_mean_xRelevant,[-10 100],fwhmVal*[1 1]);
    [relativeLoc, relativeLocPos, relativeLocNeg] = deal(X0 - peakLoc);
    relativeLocPos(relativeLoc<0) = inf;
    relativeLocNeg(relativeLoc>0) = -inf;
    
    Loc2Right = min(relativeLocPos) + peakLoc;
    Loc2Left = max(relativeLocNeg) + peakLoc;
    if nargout>1
        
        FWHM = Loc2Right-Loc2Left;
        FWHM_Y = fwhmVal;
    else    % export a stucture containing all data in case of single parameter request.
        FWHM = struct;
        FWHM.Loc2Right = min(relativeLocPos) + peakLoc;
        FWHM.Loc2Left = max(relativeLocNeg) + peakLoc;
        FWHM.FWHM = Loc2Right-Loc2Left;
        FWHM.FWHM_Y = fwhmVal;
        FWHM.peakX = peakLoc;
        FWHM.peakY = peakval;
        FWHM.AvgPhediStruct = AvgPhediStruct;
    end
end