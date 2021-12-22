function [phediOfSlopeShift, side, stretching_factor_vec] = Movie_phedi_findShiftBetweenSlopes_UseIntrpdSig(...
        signal1, signal2, PhediCoordinates,varargin)
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
    [numOOM, P] = setDefaults4function(varargin,3,0.375);
    
    %% interpolate curves
    x_sig = 1:length(signal1);
    x_evalIntrp = 1:10^-numOOM:length(signal1);
    signal1_intrpd = csaps(x_sig,signal1,P,x_evalIntrp);
    signal2_intrpd = csaps(x_sig,signal2,P,x_evalIntrp);
%     signal1_intrpd = interp1(x_sig,signal1,x_evalIntrp,'cubic');
%     signal2_intrpd = interp1(x_sig,signal2,x_evalIntrp,'cubic');
    
    %% find closest relevant peak on interpolated curve
    PhediCoordinates_interp = round((PhediCoordinates-1).*10^numOOM+1);
    
    %-- find approximate peak using otiginal signal:
    [originalPeakLoc1, side1] = Movie_phedi_find_relevant_peak(signal1, round(PhediCoordinates));
    [originalPeakLoc2, side2] = Movie_phedi_find_relevant_peak(signal1, round(PhediCoordinates));
    
    %-- find peaks of interpolated curve:
    [interped_pks1,interped_locs1] = findpeaks(signal1_intrpd,x_evalIntrp);
    [interped_pks2,interped_locs2] = findpeaks(signal2_intrpd,x_evalIntrp);
    
    %--- find closest interpolated peak to original one
    %-- signal1:
    x_diff = bsxfun(@minus,originalPeakLoc1(:),interped_locs1(:)');
    y_diff = bsxfun(@minus,signal1(originalPeakLoc1)',interped_pks1(:)');
    Eucld_dist = sqrt(x_diff.^2 + y_diff.^2);
    [~,I1] = min(Eucld_dist,[],2);
    %-- signal2:
    x_diff = bsxfun(@minus,originalPeakLoc2(:),interped_locs2(:)');
    y_diff = bsxfun(@minus,signal2(originalPeakLoc2)',interped_pks2(:)');
    Eucld_dist = sqrt(x_diff.^2 + y_diff.^2);
    [~,I2] = min(Eucld_dist,[],2);
    %-- transform to new interpolated coordinates:
    [~,interped_locs1] = findpeaks(signal1_intrpd);
    [~,interped_locs2] = findpeaks(signal2_intrpd);
    chosenPeakLoc1 = interped_locs1(I1);
    chosenPeakLoc2 = interped_locs2(I2);
    % chosenPeakLoc1 = (interped_locs1(I1)-1)*10^numOOM+1;
    % chosenPeakLoc2 = (interped_locs2(I2)-1)*10^numOOM+1;
    
    % [chosenPeakLoc1, side1] = Movie_phedi_find_relevant_peak(signal1, PhediCoordinates_interp);
    % [chosenPeakLoc2, side2] = Movie_phedi_find_relevant_peak(signal2, PhediCoordinates_interp);
    
    % [chosenPeakLoc1, side1]= Movie_phedi_find_relevant_peak(smoothed1, PhediCoordinates);
    % [chosenPeakLoc2, side2] = Movie_phedi_find_relevant_peak(smoothed2, PhediCoordinates);
    side = side2;
    side(~(side1==side2)) = 0;
    
    %% iterate upon each phedi
    
    phediOfSlopeShift = zeros(size(PhediCoordinates_interp));
    stretching_factor_vec = ones(size(PhediCoordinates_interp));
    for phedi_idx = 1:length(PhediCoordinates_interp)
        current_phedi = PhediCoordinates_interp(phedi_idx);
        currentPkLoc1 = chosenPeakLoc1(phedi_idx);
        currentPkLoc2 = chosenPeakLoc2(phedi_idx);
        %% strech second signal:
        try
            stretching_factor = signal1_intrpd(currentPkLoc1)/signal2_intrpd(currentPkLoc2);
        catch
            disp('');
        end
        %     signal2_streched = interpedSig_2*stretching_factor;
        
        %% get only section of slope:
        %-- this will include the section from the peak found and into the lower
        %part from the other side of bottomOfSlopeCoordinate.
        %     signal2location = 1:length(signal2);
        
        HorzDist2Peak = abs(current_phedi-currentPkLoc1);
        original_slope_location = current_phedi-HorzDist2Peak:current_phedi+HorzDist2Peak;
        
        %---    here each row represents a comparison for a different shift to
        %-      apply and interpolate to the second signal.
        shiftsVec = unique([-HorzDist2Peak,HorzDist2Peak,(-HorzDist2Peak:10^numOOM:HorzDist2Peak)]);
        for OO = (numOOM-1):-1:-1
            current_slope_location = current_phedi-HorzDist2Peak:current_phedi+HorzDist2Peak;
            Nshifts = length(shiftsVec);
            Npoints = length(current_slope_location);
            shiftsMat = repmat(current_slope_location,Nshifts,1)+shiftsVec(:)*ones(1,Npoints);
            %--- eliminate too-small and too large indicies
            shiftsMat(shiftsMat<1) = 1;
            shiftsMat(shiftsMat>length(signal2_intrpd)) = length(signal2_intrpd);
            %-- create matrix of signal 2 shofted in various ways
            try
                signal2_shifted = signal2_intrpd(shiftsMat);
            catch
                disp('');
            end
            slope2_mins = bsxfun(@minus,signal2_shifted,min(signal2_shifted,[],2));
            slope2_strchd = bsxfun(@rdivide,slope2_mins,max(slope2_mins,[],2));
            slope1_mins = signal1_intrpd(original_slope_location)-min(signal1_intrpd(original_slope_location));
            slope1_strchd = slope1_mins./max(slope1_mins);
            RMSofPixelSteps = rms(bsxfun(@minus,slope2_strchd,slope1_strchd),2);
            [~,minimalRMSindex] = min(RMSofPixelSteps);
            %-- update shifts vec:
            %         minimalRMSindex, shiftsVec(([-1 0 1]+minimalRMSindex))
            shiftsVec_prev = shiftsVec;
            startNew = minimalRMSindex-1; if startNew==0; startNew=1; end;
            endNew = minimalRMSindex+1; if endNew>length(shiftsVec_prev); endNew=length(shiftsVec_prev); end;
            shiftsVec = (shiftsVec_prev(startNew):10^OO:shiftsVec_prev(endNew));        % this is the shiftsVec for the next iteration
            
        end
        phediShiftIntrp = shiftsVec_prev(minimalRMSindex);
        phediShiftReal = phediShiftIntrp*(10^-numOOM);  %this is to convert the shift to units of the original pixels
        
        phediOfSlopeShift(phedi_idx) = phediShiftReal;
        stretching_factor_vec(phedi_idx) = stretching_factor;
    end
    
end