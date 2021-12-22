function [phediShift, side] = Movie_phedi_findShiftBetweenSlopes(...
        signal1, signal2, PhediCoordinates,varargin)
    %Movie_findShiftBetweenSlopes will find the shift between two slopes,
    % by measuring all the distances between points in one slope to
    % segments in the other (both ways).
    %
    % INPUTS:
    % signal1, signal2 - the rows to compare.
    % PhediCoordinates - a 1xN vector of the location of the phedis in the
    % last frame.
    %
    % OPTIONAL INPUTS:
    % numOOM - order of magnitudes to find the shift in
    % smooth - the smooth parameter.
    %
    % OUTPUTS:
    % phediShift - is the relative distance (translation in pixles) NOT
    % THE NEW POSITION
    % side - returns the sign of the slope adjacent to the specific Phedi. will
    % return '0' if there is some inconsistency in results of slope.
    %
    % update on 29-11-2018:
    %   This new version will find a common slope by taking the common
    %   pixles between the peak and the valley of the slope. Then by
    %   minimizing a functional of distances of points to segments the
    %   ideal shift will be found.
    
    
    %% set parameters
    [numOOM, slopePenaltyPower, smoothParam] = setDefaults4function(varargin,3,0,7);
%     PhediCoordinates = PhediCoordinates(1,:);
    PhediCoordinates = PhediCoordinates(:)';
    
    %% smooth signals
    signal1smth = smooth(signal1,smoothParam);
    signal2smth = smooth(signal2,smoothParam);
    
    %% get relevant regions for fitting
    expandLowBound = -1;
    expandHighBound = 0;
    [side1, sectionBounds_curr] = phedi_getSlopeBoundsAndSide(signal1smth, round(PhediCoordinates), expandLowBound, expandHighBound);
    [side2, sectionBounds_prev] = phedi_getSlopeBoundsAndSide(signal2smth, round(PhediCoordinates), expandLowBound, expandHighBound);
    
    leftBound = max([sectionBounds_curr(:,1),sectionBounds_prev(:,1)],[],2);
    rightBound = min([sectionBounds_curr(:,2),sectionBounds_prev(:,2)],[],2);
    SlopeRegions = [leftBound, rightBound];
    
    %% determine side of phedi
    side = side2;
    side(~(side1==side2)) = 0;
    
    %% iterate upon each phedi
    shiftsVec = -1:10^-numOOM:1;
    phediShift = zeros(size(PhediCoordinates));
    for phedi_idx = 1:length(PhediCoordinates)
        slope1 = signal1smth(SlopeRegions(phedi_idx,1):SlopeRegions(phedi_idx,2));
        slope2 = signal2smth(SlopeRegions(phedi_idx,1):SlopeRegions(phedi_idx,2));
        bestShift = phedi_findSlopeShiftbyMinimztn(slope1, slope2, shiftsVec,slopePenaltyPower);
        phediShift(phedi_idx) = bestShift;
    end
    
end
