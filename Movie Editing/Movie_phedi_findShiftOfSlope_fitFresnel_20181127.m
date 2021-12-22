function [newPhediLoc, side, strechFactorVert, strechFactorHorz,HorzDisplacement,VertDisplacement, polynomFit] =...
        Movie_phedi_findShiftOfSlope_fitFresnel(signal, PhediCoordinates)
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
    PhediCoordinates = PhediCoordinates(:)';
    spatialVec = 1:length(signal);
    %% find relevant section for fit
    %--- find closest peak to phedi
    [originalPeakLoc, side] = Movie_phedi_find_relevant_peak(signal, round(PhediCoordinates));
    %-- find closest minima (between phedis) to a phedi.This will define
    %-- the relevant part for a fit
    [originalMinLoc, ~] = Movie_phedi_find_relevant_peak(-smooth(signal,3), round(PhediCoordinates));
    
    %-- use peaks an mins to get the relevant sections:
    BoundHigh = originalPeakLoc(:)+side(:);
%     BoundHigh = originalPeakLoc(:);
%     BoundLow = originalMinLoc(:)+side(:);
    BoundLow = originalMinLoc(:);
    
    sectionBounds = sort([BoundHigh,BoundLow],2);
    
    %% get normalization factors for fit
    firstBoundPeakFit = originalPeakLoc(:)+side(:);
    scondBoundPeakFit = originalPeakLoc(:)-2*side(:);
    fitSections = sort([firstBoundPeakFit(:),scondBoundPeakFit(:)],2);
    
    degree = 4;
    maxX = zeros(size(PhediCoordinates));
    maxY = zeros(size(PhediCoordinates));
    polynomFit = zeros(length(PhediCoordinates),degree+1);
    for phedi_idx = 1:length(PhediCoordinates)
        xFit = spatialVec(fitSections(phedi_idx,1):fitSections(phedi_idx,2));
        yFit = signal(fitSections(phedi_idx,1):fitSections(phedi_idx,2));
        p = polyfit(xFit,yFit,degree);
        %--- find local maximum
        p_deriv = polyder(p);
        r = roots(p_deriv);
        [~,I] = min(abs(originalPeakLoc(phedi_idx) - r));
        maxX(phedi_idx) = r(I);
        maxY(phedi_idx) = polyval(p,r(I));
        
        polynomFit(phedi_idx,:) = p(:)';
    end
    
    
    %% iterate upon each phedi
    
    HorzDisplacement = zeros(size(PhediCoordinates));
    VertDisplacement = zeros(size(PhediCoordinates));
    strechFactorVert = zeros(size(PhediCoordinates));
    strechFactorHorz = zeros(size(PhediCoordinates));
    newPhediLoc =  zeros(size(PhediCoordinates));
    for phedi_idx = 1:length(PhediCoordinates)
        cropRegion = sectionBounds(phedi_idx,1):sectionBounds(phedi_idx,2);
        X_fit = cropRegion;
%         Y_fit = signal(cropRegion)./maxY(phedi_idx);
        Y_fit = signal(cropRegion);
        [XportCoeffs, fresnelComplex] = phedi_fitFresnelFunction(X_fit,Y_fit);
        strechFactorVert(phedi_idx) = XportCoeffs(1);
        strechFactorHorz(phedi_idx) = XportCoeffs(2);
        HorzDisplacement(phedi_idx) = XportCoeffs(3);
        VertDisplacement(phedi_idx) = XportCoeffs(4);
        
        %--- get new phedi location
        solveFunc = @(x) fresnelComplex(XportCoeffs,x)-interp1(spatialVec,signal,PhediCoordinates(phedi_idx));
        newPhediLoc(phedi_idx) = fzero(solveFunc,PhediCoordinates(phedi_idx));
    end
    
end