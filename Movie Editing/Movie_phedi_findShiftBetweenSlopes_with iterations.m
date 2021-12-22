function [phediOfSlopeShift, side, stretching_factor_vec] = Movie_phedi_findShiftBetweenSlopes(signal1, signal2, PhediCoordinates,varargin)
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
    isolatedSlope = signal1(original_slope_location);
    
    
    %% iterate to fit slope to second signal by minimal RMS:
    %--- find minimal RMS in pixel steps:
    pixShifts = -HorzDist2Peak:HorzDist2Peak;
    RMSofPixelSteps = zeros(size(pixShifts));
    for pixStep = 1:length(pixShifts)
        
        %--- create the bouderies for the coordinations to crop:
        tempVec = [PhediCoordinates(phedi_idx)+pixShifts(pixStep)+HorzDist2Peak,PhediCoordinates(phedi_idx)+pixShifts(pixStep)-HorzDist2Peak];
        ExtraFromLeft = nnz(min(tempVec):max(tempVec)<1);
        tempVec(tempVec<1) = 1;
        ExtraFromRight = nnz(min(tempVec):max(tempVec)>length(signal2));
        tempVec(tempVec>length(signal2)) = length(signal2);
        signal2CroppingLocations = signal2location(min(tempVec):max(tempVec));
        signal2CroppingLocations = padarray(signal2CroppingLocations, [0 ExtraFromLeft],'replicate','pre');
        signal2CroppingLocations = padarray(signal2CroppingLocations, [0 ExtraFromRight],'replicate','post');
        signal2currentSection = signal2_streched(signal2CroppingLocations);
        try
            RMSofPixelSteps(pixStep) = rms(signal2currentSection-isolatedSlope);
        catch
            disp('err');
        end
    end
    [~, minIdx] = min(RMSofPixelSteps);
    newReferencePoint = pixShifts(minIdx)+PhediCoordinates(phedi_idx);
    
    %--- find smaller shifts:
    OOM = logspace(-1,-numOOM,numOOM);
        for oomIdx = 1:numOOM
            xq = original_slope_location(1):OOM(oomIdx):original_slope_location(end);
            interpIsolatedSlope = interp1(original_slope_location,isolatedSlope,xq);
            
            pixShifts = -(OOM(oomIdx)*10):OOM(oomIdx):(OOM(oomIdx)*10);
            RMSofPixelSteps = zeros(size(pixShifts));
            for pixStep = 1:length(pixShifts)
                tempVec = [newReferencePoint+pixShifts(pixStep)+HorzDist2Peak,newReferencePoint+pixShifts(pixStep)-HorzDist2Peak];
                WantedArea4Calc = min(tempVec):OOM(oomIdx):max(tempVec);
                %         signal2CroppingLocations =
                signal2currentSection = interp1(signal2location,signal2_streched,WantedArea4Calc);
                %         signal2CroppingLocations = signal2location(min(tempWantedArea):max(tempWantedArea));
                %         signal2currentSection = signal2_streched(signal2CroppingLocations);
                RMSofPixelSteps(pixStep) = rms(signal2currentSection-interpIsolatedSlope);
            end
            [~, minIdx] = min(RMSofPixelSteps);
            newReferencePoint = pixShifts(minIdx)+newReferencePoint;
            
        end
    phediOfSlopeShift(phedi_idx) = newReferencePoint-PhediCoordinates(phedi_idx);
    stretching_factor_vec(phedi_idx) = stretching_factor;
end


end
