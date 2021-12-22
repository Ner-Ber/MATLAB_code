function [frameCount, PhediLocation, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift_updating(RowOverTime, startRow, endRow, PhediCoordinates, varargin)
%Movie_calculate_continues_shift will use the function Movie_findShiftBetweenSlopes
% to calculate the sift of slopes as function of time, taking bottomOfSlopesCoordinates
% as the bottom edge of the slope.
%
% PhediCoordinates - is a vector containing the bottoms of slopes
% in the first row (entered on the function as startRow). {Phedi in Nepali
% means "base of a hill"}
%relative2first - compare to first row or two consecutive rows
%% set defaults
[numOOM, relative2first] = setDefaults4function(varargin,3,1);

frameCount = startRow:endRow;
PhediLocation = zeros(length(frameCount),numel(PhediCoordinates));
StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
slopeIncline = cell(1,numel(PhediCoordinates));
%% run the Movie_phedi_findShiftBetweenSlopes on the Phedis given
PhediCoordinatesUpdate = PhediCoordinates(:)';
%--- find the shift between lines and the first line
if relative2first
    for rowIdx = 1:length(frameCount)
        [PhediShift, slopeIncline, stretching_factor] = Movie_phedi_findShiftBetweenSlopes(...
            RowOverTime(startRow,:),RowOverTime(frameCount(rowIdx),:), PhediCoordinates,...
            [],numOOM);
        PhediLocation(rowIdx,:) = PhediShift(:)';
        StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
    end
    %--- find the shift between lines and the first line
else
    
    
    for rowIdx = 2:length(frameCount)
        try
            [PhediShift, slopeIncline, stretching_factor] =Movie_phedi_findShiftBetweenSlopes_4updating(...
                RowOverTime(frameCount(rowIdx-1),:),RowOverTime(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
                [],numOOM);
            PhediLocation(rowIdx,:) = PhediShift(:)';
            StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
            PhediCoordinatesUpdate = PhediCoordinatesUpdate+PhediShift(:)';
        catch
            disp('');
        end
    end
    PhediLocation = cumsum(PhediLocation,1);
end
end