function [frameCount, PhediLocation, slopeIncline, StrechingFactorsMat, SlopeRegions] = Movie_phedi_calculate_continues_shift(...
        RowOverTime, startRow, endRow, PhediCoordinates, varargin)
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
    % This version was changed on the 29/11/2018.
    % Changes: The function now directs to the daughter function
    % phedi_findSlopeShiftbyMinimztn which uses a fnctional of distances,
    % weighted by the slopes of the segments fitted in order to find the
    % ideal shift between two slopes.
    %
    %
    % This version was changed on the 13/11/2018.
    % Changes: To save time, this version creates a structure of each
    % relevant row interpolation (using the function 'csaps') which is
    % inserted into the function
    % 'Movie_phedi_findShiftBetweenSlopes_UseIntrpStruct'
    %
    %
    %
    % prev. change on the 10/11/2018.
    % Changes: the rows in row over time are interpolated here, so the
    % function 'Movie_phedi_findShiftBetweenSlopes_IterMat' get
    % interpolated signals. This way it cuts in haf the number of time
    % interpolation is needed.
    % The QAD (Qiuck and dirty) method is no longer relevant for this
    % method now.
    
    %% set defaults
    [numOOM, relative2first, slopePenalty, QAD,interpMethod, P] = setDefaults4function(varargin,3,1,-1,0,'v5cubic',0.375);
    
    frameCount = startRow:endRow;
    PhediLocation = zeros(length(frameCount),numel(PhediCoordinates));
    StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
%     slopeIncline = cell(1,numel(PhediCoordinates));
    slopeIncline = [];
    SlopeRegions = [];
    %% interpolate rows:
    %-- set parameters for interpolation:
    Sig = RowOverTime(1,:);
    x_Sig = 1:length(Sig);
    Nknots = length(Sig);
    Degree = 3;
    %     x_interp = 1:10^-numOOM:length(Sig);
    %-- create meshgrods for 2D interpolation:
    %     [Xi,Ri] = meshgrid(x_Sig,1:length(frameCount));
    %     [Xq,Rq] = meshgrid(x_interp,1:length(frameCount));
    %     RowOverTime_interped = interp2(Xi,Ri,RowOverTime(frameCount,:),Xq,Rq,'spline');
    %     PhediCoordinates_interp = find(ismember(x_interp,PhediCoordinates))
    
    
    %% run the Movie_phedi_findShiftBetweenSlopes on the Phedis given
    %     PhediCoordinatesUpdate = PhediCoordinates_interp(:)';
    %     PhediCoordinatesUpdate = PhediCoordinates;
    %--- find the shift between lines and the first line
    if relative2first
        for rowIdx = 1:length(frameCount)
            [PhediShift, slopeIncline_i, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_IterMat(...
                RowOverTime_interped(1,:),RowOverTime_interped(rowIdx,:), PhediCoordinates_interp,...
                [],numOOM);
            PhediLocation(rowIdx,:) = PhediShift(:)';
            StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
            slopeIncline = cat(1,slopeIncline,slopeIncline_i);
        end
        %--- find the shift between lines and the first line
    else
        %         PhediCoordinatesUpdate = (PhediCoordinates-1)*10^numOOM+1;
        PhediCoordinatesUpdate = PhediCoordinates;
        PhediLocation(1,:) = PhediCoordinatesUpdate;
        %         StructSig1 = slmengine(x_Sig,RowOverTime(frameCount(1),:),'degree',Degree,'knots',Nknots);
        for rowIdx = 2:length(frameCount)
            %--- this uses the slm package : File Exchange - SLM - Shape
            %-- Language Modeling:
            %             StructSig2 = slmengine(x_Sig,RowOverTime(frameCount(rowIdx),:),'degree',Degree,'knots',Nknots);
            
            %             [PhediShift, slopeIncline, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_UseIntrpStruct(...
            %                 StructSig1,StructSig2, PhediCoordinatesUpdate, [],numOOM);
            
            %             x_evalPeaks = 1:10^-numOOM:length(Sig);
            %             signal1 = slmeval(x_evalPeaks,StructSig1);
            %             signal2 = slmeval(x_evalPeaks,StructSig2);
            %             [PhediShift, slopeIncline, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_UseIntrpdSig(...
            %                 signal1,signal2, PhediCoordinatesUpdate, [],numOOM);
                        [PhediShift, slopeIncline, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_UseIntrpdSig(...
                            RowOverTime(frameCount(rowIdx-1),:),RowOverTime(frameCount(rowIdx),:), PhediCoordinatesUpdate, numOOM, P);
%                         [PhediShift, slopeIncline, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_UseIntrpStruct(...
%                             StructSig1,StructSig2, PhediCoordinatesUpdate, [],numOOM);
%             [PhediShift, slopeIncline_i, SlopeRegions_i] = Movie_phedi_findShiftBetweenSlopes(...
%                 RowOverTime(frameCount(rowIdx-1),:),RowOverTime(frameCount(rowIdx),:), PhediCoordinatesUpdate, numOOM, slopePenalty);
            
%             PhediLocation(rowIdx,:) = PhediShift(:)';
%             StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
            %                 PhediCoordinatesUpdate = round(PhediCoordinatesUpdate+PhediShift(:)');
            
            PhediCoordinatesUpdate = (PhediCoordinatesUpdate-PhediShift(:));
            PhediLocation(rowIdx,:) = PhediCoordinatesUpdate;
            slopeIncline = cat(1,slopeIncline,slopeIncline_i);
            SlopeRegions = cat(3,SlopeRegions,SlopeRegions_i);
            %--- load the new interpolation for use as the previous in the
            %-- next iteration:
            %             StructSig1 = StructSig2;
        end
%         PhediLocation = cumsum(PhediLocation,1);
        %         PhediLocation = PhediLocation*10^-numOOM;
    end
    slopeIncline = mean(slopeIncline,1);
    
end