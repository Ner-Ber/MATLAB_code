function [frameCount, PhediLocation, slopeInclineMat, StrechingFactorsMat, SlopeRegions] = Movie_phedi_calculate_continues_shift_total(...
        RowOverTime_norm, startRow, endRow, PhediCoordinates, method, varargin)
    %[frameCount, PhediLocation, slopeIncline, StrechingFactorsMat, SlopeRegions] = Movie_phedi_calculate_continues_shift_total(...
    %   RowOverTime_norm, startRow, endRow, PhediCoordinates, OPTIONAL)
    %
    %Movie_phedi_calculate_continues_shift_total will use the phedi shift
    %in various ways according to user's input.
    %
    % PhediCoordinates - is a vector containing the bottoms of slopes
    % in the first row (entered on the function as startRow). {Phedi in Nepali
    % means "base of a hill"}
    %
    %   METHOD:
    %       VARARGIN:
    %   'RMS_minimization',
    %       [numOOM, slopePenalty, expandLowBound, expandHighBound, relative2first] = (4,0,-2,0,0);
    %   'csaps_intrp'
    %       [numOOM, P] = (3,0.5);
    %   'csaps_intrp_NoStrch'
    %       [numOOM, P] = (3,0.1);
    %   'iterMAT'
    %       [numOOM] = (3);
    %   'collapse_previous'
    %       [-]
    %   'collapse_prevAndFit'
    %       [NUseForPrtrb] = (100);
    %   'fit_function'
    %       none
    %
    %     
    % 29/5/2020 - added a try-catch so will not collapse if irrelevant
    % phedi found, but rather return what has already been accomplished
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
    
    %% set default method
    if isempty(method)
        method = 'RMS_minimization';
    end
    %     [numOOM, relative2first, slopePenalty, QAD,interpMethod, P, NUseForPrtrb] = setDefaults4function(varargin,3,1,-1,0,'v5cubic',0.375, 100);
    %% set universal variables
    frameCount = startRow:endRow;
    
    %% find phedi shift according to method
    if strcmpi(method,'RMS_minimization')
        %% use rms minimization
        [numOOM, slopePenalty, expandLowBound, expandHighBound, relative2first, smoothParam] = setDefaults4function(varargin,4,0,-2,0,1);
        PhediCoordinatesUpdate = PhediCoordinates;
        
        PhediLocation = zeros(length(frameCount),numel(PhediCoordinates));
        PhediLocation(1,:) = PhediCoordinatesUpdate;
        StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
        slopeInclineMat = [];
        SlopeRegions = [];
        for rowIdx = 2:length(frameCount)
            try
                if relative2first
                    [PhediShift, slopeIncline_i, SlopeRegions_i] = Movie_phedi_findShiftBetweenSlopes_RMS_relative2first(...
                        RowOverTime_norm(frameCount(1),:),RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinates, PhediCoordinatesUpdate,...
                        expandLowBound, expandHighBound, numOOM, slopePenalty, smoothParam);
                else
                    [PhediShift, slopeIncline_i, SlopeRegions_i] = Movie_phedi_findShiftBetweenSlopes_RMS(...
                        RowOverTime_norm(frameCount(rowIdx-1),:),RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
                        expandLowBound, expandHighBound, numOOM, slopePenalty, smoothParam);
                end
    %             [PhediShift, slopeIncline_i, SlopeRegions_i] = Movie_phedi_findShiftBetweenSlopes_RMS(...
    %                 RowOverTime_norm(frameCount(preIdx),:),RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
    %                 expandLowBound, expandHighBound, numOOM, slopePenalty);

                PhediCoordinatesUpdate = (PhediCoordinatesUpdate-PhediShift(:));
                PhediLocation(rowIdx,:) = PhediCoordinatesUpdate;
                slopeInclineMat = cat(1,slopeInclineMat,slopeIncline_i);
                SlopeRegions = cat(3,SlopeRegions,SlopeRegions_i);
            catch
                break
            end
            
            
        end
    elseif strcmpi(method,'iterMAT')
        %% iterate on matrix using linear interpolation
        %--- set defaults
        [numOOM] = setDefaults4function(varargin,3);
        PhediCoordinatesUpdate = PhediCoordinates(:)';
        PhediShiftMat = zeros(length(frameCount),numel(PhediCoordinates));
        PhediShiftMat(1,:) = PhediCoordinatesUpdate;
        StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
        %     slopeIncline = cell(1,numel(PhediCoordinates));
        slopeInclineMat = [];
        SlopeRegions = [];
        PhediCoordinatesUpdateMAT = [];
        StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
        for rowIdx = 2:length(frameCount)
            %             [PhediShift, slopeIncline, stretching_factor] =Movie_phedi_findShiftBetweenSlopes_IterMat(...
            %                 RowOverTime(frameCount(rowIdx-1),:),RowOverTime(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
            %                 [],numOOM);
            %             [PhediShift, slopeIncline_i, stretching_factor] =Movie_phedi_findShiftBetweenSlopes_IterMat_20181106(...
            %                 RowOverTime_norm(frameCount(rowIdx-1),:),RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
            [PhediShift, slopeIncline_i, stretching_factor] =Movie_phedi_findShiftBetweenSlopes_IterMat(...
                RowOverTime_norm(frameCount(rowIdx-1),:),RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
                [],numOOM);
            PhediShiftMat(rowIdx,:) = PhediShift(:)';
            slopeInclineMat = cat(1,slopeInclineMat,slopeIncline_i);
            StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
            PhediCoordinatesUpdate = round(PhediCoordinatesUpdate+PhediShift(:)');
            PhediCoordinatesUpdateMAT = cat(1,PhediCoordinatesUpdateMAT,PhediCoordinatesUpdate+PhediShift(:)');
        end
        PhediLocation = cumsum(PhediShiftMat,1);
        
    elseif strcmpi(method,'csaps_intrp_NoStrch')
        %% use csaps interpolation and smooth NO STRECHING
        [numOOM, P] = setDefaults4function(varargin,3,0.1);
        PhediCoordinatesUpdate = PhediCoordinates(:)';
        slopeInclineMat = [];
        PhediShiftMat = zeros(length(frameCount),numel(PhediCoordinates));
        StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
        for rowIdx = 2:length(frameCount)
            %-- set parameters for interpolation:
            %             Sig = RowOverTime_norm(1,:);
            %             x_Sig = 1:length(Sig);
            %             x_interp = 1:10^-numOOM:length(Sig);
            %             %-- create meshgrods for 2D interpolation:
            %             [Xi,Ri] = meshgrid(x_Sig,1:length(frameCount));
            %             [Xq,Rq] = meshgrid(x_interp,1:length(frameCount));
            %             RowOverTime_interped = interp2(Xi,Ri,RowOverTime_norm(frameCount,:),Xq,Rq,'spline');
            %             PhediCoordinates_interp = find(ismember(x_interp,PhediCoordinates));
            
            [PhediShift, slopeIncline_i, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_UseIntrpdSig_NoStrch(...
                RowOverTime_norm(frameCount(rowIdx-1),:),RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinatesUpdate, numOOM, P);
            PhediShiftMat(rowIdx,:) = PhediShift(:)';
            PhediCoordinatesUpdate = round(PhediCoordinatesUpdate+PhediShift(:)');
            slopeInclineMat = cat(1,slopeInclineMat,slopeIncline_i);
%             StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
            StrechingFactorsMat = [];
        end
        PhediLocation = cumsum(PhediShiftMat,1);
        SlopeRegions = [];
    
    elseif strcmpi(method,'csaps_intrp')
        %% use csaps interpolation and smooth
        [numOOM, P] = setDefaults4function(varargin,3,0.5);
        PhediCoordinatesUpdate = PhediCoordinates(:)';
        slopeInclineMat = [];
        PhediShiftMat = zeros(length(frameCount),numel(PhediCoordinates));
        StrechingFactorsMat = ones(length(frameCount),numel(PhediCoordinates));
        for rowIdx = 2:length(frameCount)
            %-- set parameters for interpolation:
            %             Sig = RowOverTime_norm(1,:);
            %             x_Sig = 1:length(Sig);
            %             x_interp = 1:10^-numOOM:length(Sig);
            %             %-- create meshgrods for 2D interpolation:
            %             [Xi,Ri] = meshgrid(x_Sig,1:length(frameCount));
            %             [Xq,Rq] = meshgrid(x_interp,1:length(frameCount));
            %             RowOverTime_interped = interp2(Xi,Ri,RowOverTime_norm(frameCount,:),Xq,Rq,'spline');
            %             PhediCoordinates_interp = find(ismember(x_interp,PhediCoordinates));
            
            [PhediShift, slopeIncline_i, stretching_factor] = Movie_phedi_findShiftBetweenSlopes_UseIntrpdSig(...
                RowOverTime_norm(frameCount(rowIdx-1),:),RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinatesUpdate, numOOM, P);
            PhediShiftMat(rowIdx,:) = PhediShift(:)';
            PhediCoordinatesUpdate = round(PhediCoordinatesUpdate+PhediShift(:)');
            slopeInclineMat = cat(1,slopeInclineMat,slopeIncline_i);
            StrechingFactorsMat(rowIdx,:) = stretching_factor(:)';
        end
        PhediLocation = cumsum(PhediShiftMat,1);
        SlopeRegions = [];
        
    elseif strcmpi(method,'fit_function')
        %% fit a fresnel or a tanh to find shift
        %-- set parameters
        
        PhediLocation = zeros(length(frameCount),numel(PhediCoordinates));
        ShiftsMat = zeros(length(frameCount)+1,length(PhediCoordinates));
        strechFactorVertMat = zeros(length(frameCount),length(PhediCoordinates));
        strechFactorHorzMat = zeros(length(frameCount),length(PhediCoordinates));
        VertDisplacementMat = zeros(length(frameCount),length(PhediCoordinates));
        slopeInclineMat = zeros(length(frameCount),length(PhediCoordinates));
        
        %--  run the Movie_phedi_findShiftOfSlope_fitFresnel on the Phedis given
        
        PhediCoordinatesUpdate = PhediCoordinates;
        rowForFindngParams = 4;
        %--- fit first rows to find parameters for the fitting function
        for rowIdx = 1:rowForFindngParams
            [newPhediLoc, side, strechFactorVert, strechFactorHorz, HorzDisplacement, VertDisplacement, polynomShift] =...
                Movie_phedi_findShiftOfSlope_fitFresnel(RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinatesUpdate);
            ShiftsMat(rowIdx+1,:) = HorzDisplacement(:)';
            strechFactorVertMat(rowIdx,:) = strechFactorVert(:)';
            strechFactorHorzMat(rowIdx,:) = strechFactorHorz(:)';
            VertDisplacementMat(rowIdx,:) = VertDisplacement(:)';
            slopeInclineMat(rowIdx,:) = side(:)';
            
            %-- update phedi location using the relative shif of the fit
            %         NewDisplacement = ShiftsMat(rowIdx+1,:)-ShiftsMat(rowIdx,:);
            %         PhediCoordinatesUpdate = newPhediLoc(:)';
            NewDisplacement = ShiftsMat(rowIdx+1,:);
            PhediCoordinatesUpdate = NewDisplacement(:)';
            PhediLocation(rowIdx,:) = PhediCoordinatesUpdate;
        end
        
        empiricStrechFactorVert = mean(strechFactorVertMat(1:rowForFindngParams,:),1);
        empiricStrechFactorVert = nan(size(empiricStrechFactorVert));
        empiricVertDisplacement = mean(VertDisplacementMat(1:rowForFindngParams,:),1);
        empiricStrechFactorHorz = mean(strechFactorHorzMat(1:rowForFindngParams,:),1);
        
        %--- use the found parameters for the fitting in the following rows
        for rowIdx = rowForFindngParams+1:length(frameCount)
            [newPhediLoc, side, strechFactorVert, strechFactorHorz, HorzDisplacement, VertDisplacement, polynomShift] =...
                Movie_phedi_findShiftOfSlope_fitFresnel(...
                RowOverTime_norm(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
                empiricStrechFactorVert, empiricStrechFactorHorz,empiricVertDisplacement);
            ShiftsMat(rowIdx+1,:) = HorzDisplacement(:)';
            strechFactorVertMat(rowIdx,:) = strechFactorVert(:)';
            strechFactorHorzMat(rowIdx,:) = strechFactorHorz(:)';
            VertDisplacementMat(rowIdx,:) = VertDisplacement(:)';
            slopeInclineMat(rowIdx,:) = side(:)';
            
            %-- update phedi location using the relative shif of the fit
            %         NewDisplacement = ShiftsMat(rowIdx+1,:)-ShiftsMat(rowIdx,:);
            %         PhediCoordinatesUpdate = newPhediLoc(:)';
            NewDisplacement = ShiftsMat(rowIdx+1,:);
            PhediCoordinatesUpdate = NewDisplacement(:)';
            PhediLocation(rowIdx,:) = PhediCoordinatesUpdate;
        end
        %--- eliminate first row of zeros:
        PhediLocation = ShiftsMat(2:end,:);
        
    elseif strcmpi(method,'collapse_previous')
        %% fit the current phedi to the collapsed previous phedis
        
        %***RIGHT A FUNCTION THAT TAKES PREVIOUS PHEDIS AND FITS THE
        %*NEW PHEDI TO THEM, MAYBE USING THE MINIMAL RMS
        
    elseif strcmpi(method,'collapse_prevAndFit')
        %% use previous phedis to build the function for the next iteration
        [NUseForPrtrb] = setDefaults4function(varargin,100);
        %--
        PhediLocation = zeros(length(frameCount)+1,numel(PhediCoordinates));
        VertStrechFactorMat = zeros(length(frameCount),length(PhediCoordinates));
        HorzSrechFactorMat = zeros(length(frameCount),length(PhediCoordinates));
        VertDisplacementMat = zeros(length(frameCount),length(PhediCoordinates));
        HorzDisplacementMat = zeros(length(frameCount),length(PhediCoordinates));
        slopeInclineMat = zeros(length(frameCount),length(PhediCoordinates));
        %---run the Movie_phedi_findShiftOfSlope_fitFresnel on the Phedis given
        PhediLocation(1,:) = PhediCoordinates(:)';
        %--- fit first rows to find parameters for the fitting function
        for rowIdx = 1:length(frameCount)
            %-- set rows to take previous data from:
            frstPertrbdFrame = max(1,rowIdx-NUseForPrtrb);
            %-- set previous phedi location:
            prevSignals = RowOverTime_norm(frameCount(frstPertrbdFrame):frameCount(rowIdx),:);
            PrevPhediData = PhediLocation(frstPertrbdFrame:rowIdx,:);
            preVertStrch = VertStrechFactorMat(frstPertrbdFrame:rowIdx,:);
            preHorzStrch = HorzSrechFactorMat(frstPertrbdFrame:rowIdx,:);
            preVertDisplcmnt = VertDisplacementMat(frstPertrbdFrame:rowIdx,:);
            preHorzDisplcmnt = HorzDisplacementMat(frstPertrbdFrame:rowIdx,:);
            %-- run the fit function
            [newPhediLoc, side, strechFactorVert, strechFactorHorz, HorzDisplacement, VertDisplacement] =...
                Movie_phedi_findShiftOfSlope_fitFunc_wPrtrbtion(prevSignals, PrevPhediData,...
                preVertStrch, preHorzStrch, preVertDisplcmnt, preHorzDisplcmnt);
            
            VertStrechFactorMat(rowIdx,:) = strechFactorVert(:)';
            HorzSrechFactorMat(rowIdx,:) = strechFactorHorz(:)';
            VertDisplacementMat(rowIdx,:) = VertDisplacement(:)';
            HorzDisplacementMat(rowIdx,:) = HorzDisplacement(:)';
            slopeInclineMat(rowIdx,:) = side(:)';
            
            %-- update phedi location using the relative shif of the fit
            %         NewDisplacement = HorzDisplacementMat(rowIdx,:);
            %         PhediCoordinatesUpdate = HorzDisplacementMat(rowIdx,:);
            PhediLocation(rowIdx+1,:) = HorzDisplacementMat(rowIdx,:);
        end
        PhediLocation = PhediLocation(2:end,:);
        
    end
    %     slopeInclineMat = mean(slopeInclineMat,1);
    
end