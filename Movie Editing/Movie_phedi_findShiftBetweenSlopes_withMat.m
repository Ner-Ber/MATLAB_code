function [phediOfSlopeShift, side, stretching_factor_vec] = Movie_phedi_findShiftBetweenSlopes_withMat(signal1, signal2, PhediCoordinates,varargin)
%Movie_findShiftBetweenSlopes will find the shift between two slopes, after
% streching tem to the same height. the height will be defined the vertical
% distance from the usr input bottomOfSlopeCoordinate to the closest peak
% (after signal is smoothed).
% bottomOfSlopeShift - is the relative distance (translation in pixles) NOT
% THE NEW POSITION
% side - returns the sign of the slope adjacent to the specific Phedi. will
% return '0' if there is some inconsistency in results of slope.
% varargin - will contain the smooth parameter


PhediCoordinates = PhediCoordinates(:)';
%% set defaults
[smoothParameter, numOOM] = setDefaults4function(varargin,6,3);
%% smooth signals
smoothed1 = smooth(signal1,smoothParameter);
smoothed2 = smooth(signal2,smoothParameter);

%% find closest relevant peak
[chosenPeakLoc1, side1]= Movie_phedi_find_relevant_peak(smoothed1, PhediCoordinates);
[chosenPeakLoc2, side2] = Movie_phedi_find_relevant_peak(smoothed2, PhediCoordinates);
side = side2;
side(~(side1==side2)) = 0;

%% iterate upon each phedi
phediOfSlopeShift = zeros(size(PhediCoordinates));
stretching_factor_vec = ones(size(PhediCoordinates));
for phedi_idx = 1:length(PhediCoordinates)
    current_phedi = PhediCoordinates(phedi_idx);
    currentPkLoc1 = chosenPeakLoc1(phedi_idx);
    currentPkLoc2 = chosenPeakLoc2(phedi_idx);
    %% strech second signal:
    stretching_factor = smoothed1(currentPkLoc1)/smoothed2(currentPkLoc2);
    signal2_streched = signal2*stretching_factor;
    
    %% get only section of slope:
    %-- this will include the section from the peak found and into the lower
    %part from the other side of bottomOfSlopeCoordinate.
    signal2location = 1:length(signal2);
    
    HorzDist2Peak = abs(current_phedi-currentPkLoc1);
    original_slope_location = current_phedi-HorzDist2Peak:current_phedi+HorzDist2Peak;
    
    
    %% iterate to fit slope to second signal by minimal RMS:
    %--- find minimal RMS in pixel steps:
    pixShifts = [-HorzDist2Peak,HorzDist2Peak];
    %--- create a vector of coordinates in which to compare to.
    %--- here each row represents a comparison for a different shift to
    %apply and interpolate to the second signal.
    tempMat = sort([PhediCoordinates(phedi_idx)+pixShifts(:)+HorzDist2Peak,PhediCoordinates(phedi_idx)+pixShifts(:)-HorzDist2Peak],2);
    C = tempMat;
    [X,Y] = meshgrid([1,10^numOOM],[1,10^numOOM]);
    [Xq,Yq] = meshgrid(1:10^numOOM,1:10^numOOM);
    Cq2 = interp2(X,Y,C,Xq,Yq);  % these are all the points in which to interpolate for comparison
    %--- eliminate too-small and too large indicies
    Cq2(Cq2<1) = 1;
    Cq2(Cq2>length(signal2)) = length(signal2);
    %--- interpolate signal 2 in order to compare
    interpedSig2 = interp1(signal2location,signal2_streched,Cq2);
    %--- interpolate signal 1 in order to compare:
    Cq1 = linspace(original_slope_location(1),original_slope_location(end),10^numOOM);
    interpedSig1 = interp1(signal2location,signal1,Cq1);
    %--- minimal RMS between the shifted and signal1:
    RMSofPixelSteps = rms(bsxfun(@minus,interpedSig2,interpedSig1),2);
    [~,minimalRMSindex] = min(RMSofPixelSteps);
    
    
    phediOfSlopeShift(phedi_idx) = Cq2(minimalRMSindex,1)-Cq1(1);
    stretching_factor_vec(phedi_idx) = stretching_factor;
end


end