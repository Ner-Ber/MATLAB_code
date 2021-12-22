function [frameCount, PhediLocation, slopeIncline, strechFactorVertMat, strechFactorHorzMat, ShiftsMat] =...
        Movie_phedi_calcShiftByFit(RowOverTime, startRow, endRow, PhediCoordinates)
    %[frameCount, PhediLocation, slopeIncline, strechFactorVertMat, strechFactorHorzMat, ShiftsMat] =...
    %   Movie_phedi_calcShiftByFit(RowOverTime, startRow, endRow, PhediCoordinates)
    %
    %Movie_phedi_calcShiftByFit will use the function Movie_phedi_findShiftOfSlope_fitFresnel
    % to calculate the sift of slopes as function of time, taking the phedi
    % coors as the bottom edge of the slope.
    %
    % PhediCoordinates - is a vector containing the bottoms of slopes
    % in the first row (entered on the function as startRow). {Phedi in Nepali
    % means "base of a hill"}
    %
    %
    % This functio bwas written on the 22-11-2018.
    % meant to replace the function
    % Movie_phedi_findShiftBetweenSlopes_UseIntrpdSig. The algorithm here
    % this fucntion will fit a complex Fresnel function to each slope and
    % calculate the shift of the function as the propagation of the phedi,
    % that is ComplexFresnel=f(x-a) 'a' is considered the shift.
    
    %% set parameters
    frameCount = startRow:endRow;
    
    PhediLocation = zeros(length(frameCount),numel(PhediCoordinates));
    ShiftsMat = zeros(length(frameCount)+1,length(PhediCoordinates));
    strechFactorVertMat = zeros(length(frameCount),length(PhediCoordinates));
    strechFactorHorzMat = zeros(length(frameCount),length(PhediCoordinates));
    VertDisplacementMat = zeros(length(frameCount),length(PhediCoordinates));
    slopeIncline = zeros(length(frameCount),length(PhediCoordinates));
    
    %% run the Movie_phedi_findShiftOfSlope_fitFresnel on the Phedis given
    
    PhediCoordinatesUpdate = PhediCoordinates;
    rowForFindngParams = 4;
    %--- fit first rows to find parameters for the fitting function
    for rowIdx = 1:rowForFindngParams
        [newPhediLoc, side, strechFactorVert, strechFactorHorz, HorzDisplacement, VertDisplacement, polynomShift] =...
            Movie_phedi_findShiftOfSlope_fitFresnel(RowOverTime(frameCount(rowIdx),:), PhediCoordinatesUpdate);
        ShiftsMat(rowIdx+1,:) = HorzDisplacement(:)';
        strechFactorVertMat(rowIdx,:) = strechFactorVert(:)';
        strechFactorHorzMat(rowIdx,:) = strechFactorHorz(:)';
        VertDisplacementMat(rowIdx,:) = VertDisplacement(:)';
        slopeIncline(rowIdx,:) = side(:)';
        
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
            RowOverTime(frameCount(rowIdx),:), PhediCoordinatesUpdate,...
            empiricStrechFactorVert, empiricStrechFactorHorz,empiricVertDisplacement);
        ShiftsMat(rowIdx+1,:) = HorzDisplacement(:)';
        strechFactorVertMat(rowIdx,:) = strechFactorVert(:)';
        strechFactorHorzMat(rowIdx,:) = strechFactorHorz(:)';
        VertDisplacementMat(rowIdx,:) = VertDisplacement(:)';
        slopeIncline(rowIdx,:) = side(:)';
        
        %-- update phedi location using the relative shif of the fit
        %         NewDisplacement = ShiftsMat(rowIdx+1,:)-ShiftsMat(rowIdx,:);
        %         PhediCoordinatesUpdate = newPhediLoc(:)';
        NewDisplacement = ShiftsMat(rowIdx+1,:);
        PhediCoordinatesUpdate = NewDisplacement(:)';
        PhediLocation(rowIdx,:) = PhediCoordinatesUpdate;
    end
    %--- eliminate first row of zeros:
    ShiftsMat = ShiftsMat(2:end,:);
end