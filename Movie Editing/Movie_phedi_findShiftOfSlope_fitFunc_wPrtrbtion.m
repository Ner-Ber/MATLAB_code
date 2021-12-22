function [newPhediLoc, side, strechFactorVert, strechFactorHorz,HorzDisplacement,VertDisplacement] =...
        Movie_phedi_findShiftOfSlope_fitFunc_wPrtrbtion(prevSignals, PrevPhediData,...
            preVertStrch, preHorzStrch, preVertDisplcmnt, preHorzDisplcmnt, varargin)
    %[newPhediLoc, side, strechFactorVert, strechFactorHorz,HorzDisplacement,VertDisplacement, polynomFit] =...
    %  Movie_phedi_findShiftOfSlope_fitFresnel(signal, PhediCoordinates, varargin)
    %
    % OPTIONAL constrains: strechFactorVertIN, strechFactorHorzIN, HorzDisplacementIN, VertDisplacementIN
    
    
    %Movie_findShiftBetweenSlopes will find the shift between two slopes, after
    % streching tem to the same height. the height will be defined the vertical
    % distance from the usr input bottomOfSlopeCoordinate to the closest peak
    % (after signal is smoothed).
    % bottomOfSlopeShift - is the relative distance (translation in pixles) NOT
    % THE NEW POSITION
    % side - returns the sign of the slope adjacent to the specific Phedi. will
    % return '0' if there is some inconsistency in results of slope.
    % varargin - will contain the smooth parameter
    
    
    %% set parameters
    %-- separate previous and current signal and fitting parameters
    PhediCoordinates = PrevPhediData(end,:);
    signal = prevSignals(end,:);
    PertrbdSignals = prevSignals(1:end-1,:);
    preVertStrch = preVertStrch(1:end-1,:);
    preHorzStrch = preHorzStrch(1:end-1,:);
    preVertDisplcmnt = preVertDisplcmnt(1:end-1,:);
    preHorzDisplcmnt = preHorzDisplcmnt(1:end-1,:);
    
    [strechFactorVertIN, strechFactorHorzIN, HorzDisplacementIN, VertDisplacementIN] = setDefaults4function(varargin,...
        nan(size(PhediCoordinates)),nan(size(PhediCoordinates)),nan(size(PhediCoordinates)),nan(size(PhediCoordinates)));
    
    spatialVec = 1:length(signal);
    
    %% find relevant section for fit
%     %--- find closest peak to phedi
%     [originalPeakLoc, side] = Movie_phedi_find_relevant_peak(signal, round(PhediCoordinates));
%     %-- find closest minima (between phedis) to a phedi.This will define
%     %-- the relevant part for a fit
%     [originalMinLoc, ~] = Movie_phedi_find_relevant_peak(-smooth(signal,3), round(PhediCoordinates));
%     
%     %-- use peaks an mins to get the relevant sections:
% %     BoundHigh = originalPeakLoc(:)+side(:);
%     BoundHigh = originalPeakLoc(:);
%     BoundLow = originalMinLoc(:)+side(:);
% %     BoundLow = originalMinLoc(:);
%     
%     sectionBounds = sort([BoundHigh,BoundLow],2);
    
    expandLowBound = 0;
    expandHighBound = 0;
    [side, sectionBounds] = phedi_getSlopeBoundsAndSide(signal, PhediCoordinates, expandLowBound, expandHighBound);
    
    
    %% iterate upon each phedi
    
    HorzDisplacement = zeros(size(PhediCoordinates));
    VertDisplacement = zeros(size(PhediCoordinates));
    strechFactorVert = zeros(size(PhediCoordinates));
    strechFactorHorz = zeros(size(PhediCoordinates));
    newPhediLoc =  zeros(size(PhediCoordinates));
    for phedi_idx = 1:length(PhediCoordinates)
        %-- first fit: find approximate fit of current phedi:
        cropRegion = sectionBounds(phedi_idx,1):sectionBounds(phedi_idx,2);
        X_fit = cropRegion;
        Y_fit = signal(cropRegion);
        [XportCoeffs1stFit, ~] = phedi_fitFresnelFunction(X_fit,Y_fit,...
            strechFactorVertIN(phedi_idx), strechFactorHorzIN(phedi_idx), HorzDisplacementIN(phedi_idx) ,VertDisplacementIN(phedi_idx));
%         [XportCoeffs1stFit, ~] = phedi_fitTanhFunction(X_fit,Y_fit,...
%             strechFactorVertIN(phedi_idx), strechFactorHorzIN(phedi_idx), HorzDisplacementIN(phedi_idx) ,VertDisplacementIN(phedi_idx));
        
        %-- collapse all data together using fitting parameters:
        cumulatedX = X_fit(:);
        cumulatedY = Y_fit(:);
        if ~isempty(PertrbdSignals)
            for prvRowIdx = 1:size(PertrbdSignals,1)
                [side, sectionBoundsPRV] = phedi_getSlopeBoundsAndSide(...
                    PertrbdSignals(prvRowIdx,:),PrevPhediData(prvRowIdx,phedi_idx), expandLowBound, expandHighBound);
                thisX = sectionBoundsPRV(1):sectionBoundsPRV(2);
                thisY = PertrbdSignals(prvRowIdx,thisX);
                
                %-- shift to zero:
                X_prev2neut = thisX- preHorzDisplcmnt(prvRowIdx,phedi_idx);
%                 X_prev2neut = (thisX- preHorzDisplcmnt(prvRowIdx,phedi_idx)).*preHorzStrch(prvRowIdx,phedi_idx);
%                 Y_prev2neut = (thisY-preVertDisplcmnt(prvRowIdx,phedi_idx))./preVertStrch(prvRowIdx,phedi_idx);
                
                %-- shift to caollapse:
                thisX_fixed = X_prev2neut + XportCoeffs1stFit(3);
                thisY_fixed = thisY;
%                 thisX_fixed = X_prev2neut/XportCoeffs1stFit(2) + XportCoeffs1stFit(3);
%                 thisY_fixed = Y_prev2neut*XportCoeffs1stFit(1) + XportCoeffs1stFit(4);                
%                 thisX_fixed = ((thisX-XportCoeffs1stFit(3))*XportCoeffs1stFit(2))/preHorzStrch(prvRowIdx,phedi_idx)+preHorzDisplcmnt(prvRowIdx,phedi_idx);
%                 thisY_fixed = ((thisY-XportCoeffs1stFit(4))/XportCoeffs1stFit(1))*preVertStrch(prvRowIdx,phedi_idx)+preVertDisplcmnt(prvRowIdx,phedi_idx);

                cumulatedX = cat(1,cumulatedX,thisX_fixed(:));
                cumulatedY = cat(1,cumulatedY,thisY_fixed(:));
            end
        end
        
        %-- fit to cumulated data
        [XportCoeffsFinalFit, ~] = phedi_fitFresnelFunction(cumulatedX,cumulatedY,...
            strechFactorVertIN(phedi_idx), strechFactorHorzIN(phedi_idx), HorzDisplacementIN(phedi_idx) ,VertDisplacementIN(phedi_idx));
        
        strechFactorVert(phedi_idx) = XportCoeffsFinalFit(1);
        strechFactorHorz(phedi_idx) = XportCoeffsFinalFit(2);
        HorzDisplacement(phedi_idx) = XportCoeffsFinalFit(3);
        VertDisplacement(phedi_idx) = XportCoeffsFinalFit(4);
        
    end
    
end