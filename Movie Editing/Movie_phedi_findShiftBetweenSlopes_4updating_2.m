function [phediOfSlopeShift, side, stretching_factor_vec] = Movie_phedi_findShiftBetweenSlopes_4updating_2(signal1, signal2, PhediCoordinates,varargin)
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
useSmoothed = 1;
%% smooth signals
smoothed1 = smooth(signal1,smoothParameter);
smoothed2 = smooth(signal2,smoothParameter);

%% interpolate signals by OOM
original_pixels = 1:length(smoothed1);
interpolate_pixel = (10^(-numOOM)):(10^(-numOOM)):length(smoothed1);
PhediCoordinatesInterp = round(PhediCoordinates*10^numOOM);
smoothedInterp1 = interp1(original_pixels,smoothed1,interpolate_pixel);
smoothedInterp2 = interp1(original_pixels,smoothed2,interpolate_pixel);
signalInterp1 = interp1(original_pixels,signal1,interpolate_pixel);
signalInterp2 = interp1(original_pixels,signal2,interpolate_pixel);
%#$#$#$#$#$ NOTICE! smoothedInterp1/2 give NaNs where there is
%extrapolation needed!!! $#$#$#$#$#%

%% find closest relevant peak
[chosenPeakLoc1, side1]= Movie_phedi_find_relevant_peak(smoothedInterp1, PhediCoordinatesInterp);
[chosenPeakLoc2, side2] = Movie_phedi_find_relevant_peak(smoothedInterp2, PhediCoordinatesInterp);
side = side2;
side(~(side1==side2)) = 0;

%% iterate upon each phedi
phediOfSlopeShiftNewGrid = zeros(size(PhediCoordinatesInterp));
stretching_factor_vec = ones(size(PhediCoordinatesInterp));
for phedi_idx = 1:length(PhediCoordinatesInterp)
    current_phedi = PhediCoordinatesInterp(phedi_idx);
    currentPkLoc1 = chosenPeakLoc1(phedi_idx);
    currentPkLoc2 = chosenPeakLoc2(phedi_idx);
    %% strech second signal:
    if useSmoothed
        stretching_factor = smoothedInterp1(currentPkLoc1)/smoothedInterp2(currentPkLoc2);
        signal2_streched = smoothedInterp2*stretching_factor;
    else
        stretching_factor = signalInterp1(currentPkLoc1)/signalInterp2(currentPkLoc2);
        signal2_streched = signalInterp2*stretching_factor;
    end
    
    
    %% use coarse shift before applying fit algorithm
    
    
    %% get only section of slope:
    %-- this will include the section from the peak found and into the lower
    %part from the other side of bottomOfSlopeCoordinate.
    signal2location = 1:length(signal2_streched);
    
    HorzDist2Peak = abs(current_phedi-currentPkLoc1);
    original_slope_location = current_phedi-HorzDist2Peak:current_phedi+HorzDist2Peak;
    if useSmoothed
        isolatedSlope = smoothedInterp1(original_slope_location);
    else
        isolatedSlope = signalInterp1(original_slope_location);
    end
    
    
    %% iterate to fit slope to second signal by minimal RMS:
    %--- find minimal RMS in pixel steps:
    newReferencePoint = PhediCoordinatesInterp(phedi_idx);
    
    OOM = logspace(numOOM,0,numOOM+1);
    OOM_LR = fliplr(OOM);
    for oomIdx = 1:length(OOM);
        shiftAmplitude = HorzDist2Peak./OOM_LR(oomIdx);
        pixShifts = round(-shiftAmplitude:OOM(oomIdx):shiftAmplitude);
        RMSofPixelSteps = zeros(size(pixShifts));
        
        for pixStep = 1:length(pixShifts)
            
            %--- create the bouderies for the coordinations to crop:
            tempVec = [newReferencePoint+pixShifts(pixStep)+HorzDist2Peak,newReferencePoint+pixShifts(pixStep)-HorzDist2Peak];
            ExtraFromLeft = nnz(min(tempVec):max(tempVec)<1);
            tempVec(tempVec<1) = 1;
            ExtraFromRight = nnz(min(tempVec):max(tempVec)>length(signalInterp2));
            tempVec(tempVec>length(signalInterp2)) = length(signalInterp2);
            if mod(min(tempVec),1)~=0 || mod(max(tempVec),1)~=0
                disp('');
            end
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
        newReferencePoint = pixShifts(minIdx)+newReferencePoint;
    end
    
    %--- find smaller shifts:
    %     OOM = logspace(-1,-numOOM,numOOM);
    %     for oomIdx = 1:numOOM
    %         xq = original_slope_location(1):OOM(oomIdx):original_slope_location(end);
    %         interpIsolatedSlope = interp1(original_slope_location,isolatedSlope,xq);
    %
    %         pixShifts = -(OOM(oomIdx)*10):OOM(oomIdx):(OOM(oomIdx)*10);
    %         RMSofPixelSteps = zeros(size(pixShifts));
    %         for pixStep = 1:length(pixShifts)
    %             tempVec = [newReferencePoint+pixShifts(pixStep)+HorzDist2Peak,newReferencePoint+pixShifts(pixStep)-HorzDist2Peak];
    %             WantedArea4Calc = min(tempVec):OOM(oomIdx):max(tempVec);
    %             %         signal2CroppingLocations =
    %             signal2currentSection = interp1(signal2location,signal2_streched,WantedArea4Calc);
    %             %         signal2CroppingLocations = signal2location(min(tempWantedArea):max(tempWantedArea));
    %             %         signal2currentSection = signal2_streched(signal2CroppingLocations);
    %             RMSofPixelSteps(pixStep) = rms(signal2currentSection-interpIsolatedSlope);
    %         end
    %         [~, minIdx] = min(RMSofPixelSteps);
    %         newReferencePoint = pixShifts(minIdx)+newReferencePoint;
    %
    %     end
    phediOfSlopeShiftNewGrid(phedi_idx) = newReferencePoint-PhediCoordinatesInterp(phedi_idx);
    stretching_factor_vec(phedi_idx) = stretching_factor;
end
phediOfSlopeShift = phediOfSlopeShiftNewGrid*(10^(-numOOM));


end