function [frameCount, PhediLocation, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift_13112018(RowOverTime, startRow, endRow, PhediCoordinates, varargin)
    %Movie_calculate_continues_shift will use the function Movie_findShiftBetweenSlopes
    % to calculate the sift of slopes as function of time, taking bottomOfSlopesCoordinates
    % as the bottom edge of the slope.
    %
    % PhediCoordinates - is a vector containing the bottoms of slopes
    % in the first row (entered on the function as startRow). {Phedi in Nepali
    % means "base of a hill"}
    %relative2first - compare to first row or two consecutive rows
    %QAD - stands for 'quick and dirty' uses a fast but inacurrate algorithm.
    %
    %
    % This version was changed on the 10/11/2018. 
    % Changes: the rows in row over time are interpolated here, so the
    % function 'Movie_phedi_findShiftBetweenSlopes_IterMat' get
    % interpolated signals. This way it cuts in haf the number of time
    % interpolation is needed. 
    % The QAD (Qiuck and dirty) method is no longer relevant for this
    % method now. 
    
    %% set defaults
    [numOOM, relative2first, QAD,interpMethod] = setDefaults4function(varargin,3,1,0,'v5cubic');
    
    frameCount = startRow:endRow;
    PhediLocation = zeros(length(frameCount),numel(PhediCoordinates));
    StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
    slopeIncline = cell(1,numel(PhediCoordinates));
    
    %% interpolate rows:
    Sig = RowOverTime(1,:);
    x_measure = 1:length(Sig);
    x_interp = 1:10^-numOOM:length(Sig);
    %-- create meshgrods for 2D interpolation:
    [Xi,Ri] = meshgrid(x_measure,1:length(frameCount));
    [Xq,Rq] = meshgrid(x_interp,1:length(frameCount));
    RowOverTime_interped = interp2(Xi,Ri,RowOverTime(frameCount,:),Xq,Rq,'spline');
    PhediCoordinates_interp = find(ismember(x_interp,PhediCoordinates));
    
    
    %% run the Movie_phedi_findShiftBetweenSlopes on the Phedis given
    PhediCoordinatesUpdate = PhediCoordinates_interp(:)';
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
                    RowOverTime_interped(1,:),RowOverTime_interped(rowIdx,:), PhediCoordinates_interp,...
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
                    RowOverTime_interped(rowIdx-1,:),RowOverTime_interped(rowIdx,:), PhediCoordinatesUpdate,...
                    [],numOOM);
                PhediLocation(rowIdx,:) = PhediShift(:)';
                StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
%                 PhediCoordinatesUpdate = round(PhediCoordinatesUpdate+PhediShift(:)');
                PhediCoordinatesUpdate = (PhediCoordinatesUpdate+PhediShift(:)'*(10^numOOM));
            end
        end
        PhediLocation = cumsum(PhediLocation,1);
    end
end