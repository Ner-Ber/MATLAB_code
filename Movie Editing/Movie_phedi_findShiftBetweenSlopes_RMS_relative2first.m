function [phediShift, side, SlopeRegions2] = Movie_phedi_findShiftBetweenSlopes_RMS_relative2first(...
        signal1, signal2, PhediCoordinatesInitial, PhediCoordinates,varargin)
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
    % expandLowBound, expandHighBound - number of points relative to
    %   minimum and maximum near phedi.
    % slopePenaltyPower - default 0 (no penalty);
    % numOOM - order of magnitudes to find the shift in
    % smoothParam - the smooth parameter.
    %
    % OUTPUTS:
    % phediShift - is the relative distance (translation in pixles) NOT
    % THE NEW POSITION
    % side - returns the sign of the slope adjacent to the specific Phedi. will
    % return '0' if there is some inconsistency in results of slope.
    %
    % update on 02-12-2018:
    %    Use interpolation to compute distances between the two slopes.
    %
    % update on 29-11-2018:
    %   This new version will find a common slope by taking the common
    %   pixles between the peak and the valley of the slope. Then by
    %   minimizing a functional of distances of points to segments the
    %   ideal shift will be found.
    
    
    %% set parameters
    [expandLowBound, expandHighBound, numOOM, slopePenaltyPower, smoothParam] = setDefaults4function(varargin,-2,0,3,0,1);
    %     PhediCoordinates = PhediCoordinates(1,:);
    PhediCoordinates = PhediCoordinates(:)';
    PhediCoordinatesInitial = PhediCoordinatesInitial(:)';
    
    %% smooth signals
    signal1smth = smooth(signal1,smoothParam);
    signal2smth = smooth(signal2,smoothParam);
    
    %% get relevant regions for fitting
    [side1, sectionBounds1, bH1, ~] = phedi_getSlopeBoundsAndSide(signal1smth, round(PhediCoordinates), expandLowBound, expandHighBound);
    [side2, sectionBounds2, bH2, ~] = phedi_getSlopeBoundsAndSide(signal2smth, round(PhediCoordinatesInitial), expandLowBound, expandHighBound);
    minL = min(abs([diff(sectionBounds1,1,2), diff(sectionBounds2,1,2)]));
    SlopeRegions1 = sort([bH1(:), bH1(:)-side1(:).*minL(:)],2);
    SlopeRegions2 = sort([bH2(:), bH2(:)-side2(:).*minL(:)],2);
    
%     leftBound = max([sectionBounds1(:,1),sectionBounds2(:,1)],[],2);
%     rightBound = min([sectionBounds1(:,2),sectionBounds2(:,2)],[],2);
%     SlopeRegions = [leftBound, rightBound];
    
    %% determine side of phedi
    side = side2;
    side(~(side1==side2)) = 0;
    
    %% iterate upon each phedi
    %     shiftsVec = -1:10^-numOOM:1;
    phediShift = zeros(size(PhediCoordinates));
    for phedi_idx = 1:length(PhediCoordinates)
        slope1 = signal1smth(SlopeRegions1(phedi_idx,1):SlopeRegions1(phedi_idx,2));
        slope2 = signal2smth(SlopeRegions2(phedi_idx,1):SlopeRegions2(phedi_idx,2));
        
%         bestShiftIntrpd = phedi_findIntrpdSlopeShiftbyMinimztn(slope1, slope2, numOOM,slopePenaltyPower);
        bestShiftIntrpd = phedi_findIntrpdSlopeShiftbyMinimztnRMS(slope1, slope2, numOOM,slopePenaltyPower);
        
%         bestShiftIntrpd = phedi_findIntrpdSlopeShiftbyMinimztnRMS_noWeight(slope1, slope2, numOOM,slopePenaltyPower);
%         bestShiftIntrpd = phedi_findIntrpdSlopeShiftbyMinimztn_orgnl2intrpd(slope1, slope2, numOOM,slopePenaltyPower);
        bestShift = bestShiftIntrpd*10^-numOOM;
%         bestShift = phedi_findSlopeShiftbyMinimztn(slope1, slope2, shiftsVec,slopePenaltyPower);
        phediShift(phedi_idx) = bestShift-(SlopeRegions2(phedi_idx,1)-SlopeRegions1(phedi_idx,1));
    end
    
end
