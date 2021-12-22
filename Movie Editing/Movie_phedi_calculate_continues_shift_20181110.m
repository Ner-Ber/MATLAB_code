function [frameCount, PhediLocation, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift_20181110(RowOverTime, startRow, endRow, PhediCoordinates, varargin)
    %Movie_calculate_continues_shift will use the function Movie_findShiftBetweenSlopes
    % to calculate the sift of slopes as function of time, taking bottomOfSlopesCoordinates
    % as the bottom edge of the slope.
    %
    % PhediCoordinates - is a vector containing the bottoms of slopes
    % in the first row (entered on the function as startRow). {Phedi in Nepali
    % means "base of a hill"}
    %relative2first - compare to first row or two consecutive rows
    %QAD - stands for 'quick and dirty' uses a fast but inacurrate algorithm.
    
    %% set defaults
    [numOOM, relative2first, QAD] = setDefaults4function(varargin,3,1,0);
    
    frameCount = startRow:endRow;
    PhediLocation = zeros(length(frameCount),numel(PhediCoordinates));
    StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
    slopeIncline = cell(1,numel(PhediCoordinates));
    %% run the Movie_phedi_findShiftBetweenSlopes on the Phedis given
    PhediCoordinatesUpdate = PhediCoordinates(:)';
    %--- find the shift between lines and the first line
    if relative2first
        if QAD
            for rowIdx = 1:length(frameCount)
                [PhediShift, slopeIncline, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_QuickAndDirty(...
                    RowOverTime(startRow,:),RowOverTime(frameCount(rowIdx),:), PhediCoordinates,...
                    [],numOOM);
                PhediLocation(rowIdx,:) = PhediShift(:)';
                StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
            end
        else
            for rowIdx = 1:length(frameCount)
                [PhediShift, slopeIncline, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_IterMat(...
                    RowOverTime(startRow,:),RowOverTime(frameCount(rowIdx),:), PhediCoordinates,...
                    [],numOOM);
                PhediLocation(rowIdx,:) = PhediShift(:)';
                StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
            end
        end
        %--- find the shift between lines and the first line
    else
        if QAD
            [PhediShift, slopeIncline, stretching_factor] =Movie_phedi_findShiftBetweenSlopes_QuickAndDirty(...
                RowOverTime(frameCount(rowIdx-1),:),RowOverTime(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
                [],numOOM);
            PhediLocation(rowIdx,:) = PhediShift(:)';
            StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
            PhediCoordinatesUpdate = round(PhediCoordinatesUpdate+PhediShift(:)');
        else
            for rowIdx = 2:length(frameCount)
                [PhediShift, slopeIncline, stretching_factor] =Movie_phedi_findShiftBetweenSlopes_IterMat(...
                    RowOverTime(frameCount(rowIdx-1),:),RowOverTime(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
                    [],numOOM);
                PhediLocation(rowIdx,:) = PhediShift(:)';
                StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
                PhediCoordinatesUpdate = round(PhediCoordinatesUpdate+PhediShift(:)');
            end
        end
        PhediLocation = cumsum(PhediLocation,1);
    end
end