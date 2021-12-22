function [side, sectionBounds,BoundHigh,BoundLow] = phedi_getSlopeBoundsAndSide(signal, PhediCoordinates, varargin)
%     [side, sectionBounds] = phedi_getSlopeBoundsAndSide(signal, PhediCoordinates, varargin)
%
%OPTIONAL:  expandLowBound, expandHighBound (must be integers)
% 
% phedi_getSlopeBoundsAndSide appears in 'Movie_phedi_findShiftOfSlope_XXX'
% functions.
    
    
    [expandLowBound, expandHighBound] = setDefaults4function(varargin,1,0);
     %--- find closest peak to phedi
    [originalPeakLoc, side] = Movie_phedi_find_relevant_peak(signal, round(PhediCoordinates));
    %-- find closest minima (between phedis) to a phedi.This will define
    %-- the relevant part for a fit
    [originalMinLoc, ~] = Movie_phedi_find_relevant_peak(-smooth(signal,3), round(PhediCoordinates));
    
    %-- use peaks an mins to get the relevant sections:
    BoundHigh = originalPeakLoc(:)+expandHighBound*side(:);
    BoundLow = originalMinLoc(:)-expandLowBound*side(:);
    
    sectionBounds = sort([BoundHigh,BoundLow],2);
    
end