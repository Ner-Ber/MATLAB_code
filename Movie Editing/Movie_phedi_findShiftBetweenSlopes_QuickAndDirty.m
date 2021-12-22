function [phediOfSlopeShift, side, stretching_factor_vec] = Movie_phedi_findShiftBetweenSlopes_QuickAndDirty(signal1, signal2, PhediCoordinates,varargin)
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
    FixedPhedi = current_phedi;
    %% iterate upon orders of magnitufe
    %--- order indicates jumps by 10^-oomIdx of a pixle
    for currentIdx = 0:numOOM
        %% iterate to fit slope to second signal by minimal RMS:
        %--- find minimal RMS in pixel steps:
        if currentIdx==0
            shiftsVec = [-HorzDist2Peak:HorzDist2Peak];
        else
            shiftsVec = [-10:10]*10^-currentIdx;
        end
        %--- here each row represents a comparison for a different shift to
        %apply and interpolate to the second signal.
        current_slope_location = FixedPhedi-HorzDist2Peak:FixedPhedi+HorzDist2Peak;
        Nshifts = length(shiftsVec);
        Npoints = length(current_slope_location);
        shiftsMat = repmat(current_slope_location,Nshifts,1)+shiftsVec(:)*ones(1,Npoints);
        %--- eliminate too-small and too large indicies
        shiftsMat(shiftsMat<1) = 1;
        shiftsMat(shiftsMat>length(signal2)) = length(signal2);
        %--- interpolate signal 2 in order to compare
        interpedSig2 = interp1(signal2location,signal2_streched,shiftsMat);
        %--- interpolate signal 1 in order to compare:
        interpedSig1 = signal1(original_slope_location);
        %--- minimal RMS between the shifted and signal1:
        RMSofPixelSteps = rms(bsxfun(@minus,interpedSig2,interpedSig1),2);
        [~,minimalRMSindex] = min(RMSofPixelSteps);
        iterationShift = shiftsVec(minimalRMSindex);
        FixedPhedi = FixedPhedi + iterationShift;
    end
    phediOfSlopeShift(phedi_idx) = FixedPhedi - current_phedi;
    stretching_factor_vec(phedi_idx) = stretching_factor;
end


end