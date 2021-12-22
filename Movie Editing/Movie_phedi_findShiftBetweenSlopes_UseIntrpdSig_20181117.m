function [phediOfSlopeShift, side, stretching_factor_vec] = Movie_phedi_findShiftBetweenSlopes_UseIntrpdSig(signal1, signal2, PhediCoordinates,varargin)
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
[smoothParameter, numOOM, interpMethod] = setDefaults4function(varargin,6,3, 'v5cubic');
% Nconvrg = 3;
%% smooth signals
% smoothed1 = smooth(signal1,smoothParameter);
% smoothed2 = smooth(signal2,smoothParameter);

%% interpolate signals
% x_measure = 1:length(signal1);
% x_interp = 1:10^-numOOM:length(signal1);
% interpedSig_1 = interp1(x_measure,signal1,x_interp,interpMethod);
% interpedSig_2 = interp1(x_measure,signal2,x_interp,interpMethod);

%% find closest relevant peak on interpolated curve
% PhediCoordinates_interp = find(ismember(x_interp,PhediCoordinates));
% x_evalPeaks = SigStruct1.x(1):10^-numOOM:SigStruct1.x(end);
% signal1 = slmeval(x_evalPeaks,SigStruct1);
% signal2 = slmeval(x_evalPeaks,SigStruct2);

% x_evalPeaks = 1:10^-numOOM:length(signal1);
% signal1 = slmeval(x_evalPeaks,SigStruct1);
% signal2 = slmeval(x_evalPeaks,SigStruct2);
% [~,PhediCoordinates] = min(abs(bsxfun(@minus,x_evalPeaks(:),PhediCoordinates(:)')),[],1);
[chosenPeakLoc1, side1] = Movie_phedi_find_relevant_peak(signal1, PhediCoordinates);
[chosenPeakLoc2, side2] = Movie_phedi_find_relevant_peak(signal2, PhediCoordinates);

% [chosenPeakLoc1, side1]= Movie_phedi_find_relevant_peak(smoothed1, PhediCoordinates);
% [chosenPeakLoc2, side2] = Movie_phedi_find_relevant_peak(smoothed2, PhediCoordinates);
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
    stretching_factor = signal1(currentPkLoc1)/signal2(currentPkLoc2);
%     signal2_streched = interpedSig_2*stretching_factor;
    
    %% get only section of slope:
    %-- this will include the section from the peak found and into the lower
    %part from the other side of bottomOfSlopeCoordinate.
%     signal2location = 1:length(signal2);
    
    HorzDist2Peak = abs(current_phedi-currentPkLoc1);
    original_slope_location = current_phedi-HorzDist2Peak:current_phedi+HorzDist2Peak;
    
    %---    here each row represents a comparison for a different shift to
    %-      apply and interpolate to the second signal.
    shiftsVec = (-HorzDist2Peak:10^numOOM:HorzDist2Peak);
    for OO = (numOOM-1):-1:-1
        current_slope_location = current_phedi-HorzDist2Peak:current_phedi+HorzDist2Peak;
        Nshifts = length(shiftsVec);
        Npoints = length(current_slope_location);
        shiftsMat = repmat(current_slope_location,Nshifts,1)+shiftsVec(:)*ones(1,Npoints);
        %--- eliminate too-small and too large indicies
        shiftsMat(shiftsMat<1) = 1;
        shiftsMat(shiftsMat>length(signal2)) = length(signal2);
        %-- create matrix of signal 2 shofted in various ways
        signal2_shifted = signal2(shiftsMat);
        slope2_mins = bsxfun(@minus,signal2_shifted,min(signal2_shifted,[],2));
        slope2_strchd = bsxfun(@rdivide,slope2_mins,max(slope2_mins,[],2));
        slope1_mins = signal1(original_slope_location)-min(signal1(original_slope_location));
        slope1_strchd = slope1_mins./max(slope1_mins);
        RMSofPixelSteps = rms(bsxfun(@minus,slope2_strchd,slope1_strchd),2);
        [~,minimalRMSindex] = min(RMSofPixelSteps);
        %-- update shifts vec:
%         minimalRMSindex, shiftsVec(([-1 0 1]+minimalRMSindex))
        shiftsVec_prev = shiftsVec;
        try
        shiftsVec = (shiftsVec(minimalRMSindex-1):10^OO:shiftsVec(minimalRMSindex+1));
        catch
           disp(); 
        end
    end
    phediShiftIntrp = shiftsVec_prev(minimalRMSindex);
%     phediShiftReal = phediShiftIntrp*(10^-numOOM);  %this is to convert the shift to units of the original pixels
    phediShiftReal = phediShiftIntrp;
  
    phediOfSlopeShift(phedi_idx) = phediShiftReal;
    stretching_factor_vec(phedi_idx) = stretching_factor;
end


end