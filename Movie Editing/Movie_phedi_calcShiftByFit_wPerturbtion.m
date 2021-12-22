function [frameCount, PhediLocation, slopeIncline, VertStrechFactorMat, HorzSrechFactorMat, HorzDisplacementMat] =...
        Movie_phedi_calcShiftByFit_wPerturbtion(RowOverTime, startRow, endRow, PhediCoordinates, varargin)
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
    [NUseForPrtrb] = setDefaults4function(varargin,100);
    frameCount = startRow:endRow;
    
    PhediLocation = zeros(length(frameCount)+1,numel(PhediCoordinates));
    VertStrechFactorMat = zeros(length(frameCount),length(PhediCoordinates));
    HorzSrechFactorMat = zeros(length(frameCount),length(PhediCoordinates));
    VertDisplacementMat = zeros(length(frameCount),length(PhediCoordinates));
    HorzDisplacementMat = zeros(length(frameCount),length(PhediCoordinates));
    slopeIncline = zeros(length(frameCount),length(PhediCoordinates));
    
    %% run the Movie_phedi_findShiftOfSlope_fitFresnel on the Phedis given
    
    PhediLocation(1,:) = PhediCoordinates(:)';
    %--- fit first rows to find parameters for the fitting function
    for rowIdx = 1:length(frameCount)
        %-- set rows to take previous data from:
        frstPertrbdFrame = max(1,rowIdx-NUseForPrtrb);
        %-- set previous phedi location:
        prevSignals = RowOverTime(frameCount(frstPertrbdFrame):frameCount(rowIdx),:);
        PrevPhediData = PhediLocation(frstPertrbdFrame:rowIdx,:);
        preVertStrch = VertStrechFactorMat(frstPertrbdFrame:rowIdx,:);
        preHorzStrch = HorzSrechFactorMat(frstPertrbdFrame:rowIdx,:);
        preVertDisplcmnt = VertDisplacementMat(frstPertrbdFrame:rowIdx,:);
        preHorzDisplcmnt = HorzDisplacementMat(frstPertrbdFrame:rowIdx,:);
        
        try
        [newPhediLoc, side, strechFactorVert, strechFactorHorz, HorzDisplacement, VertDisplacement] =...
            Movie_phedi_findShiftOfSlope_fitFunc_wPrtrbtion(prevSignals, PrevPhediData,...
            preVertStrch, preHorzStrch, preVertDisplcmnt, preHorzDisplcmnt);
        catch
            disp('');
        end
        
        VertStrechFactorMat(rowIdx,:) = strechFactorVert(:)';
        HorzSrechFactorMat(rowIdx,:) = strechFactorHorz(:)';
        VertDisplacementMat(rowIdx,:) = VertDisplacement(:)';
        HorzDisplacementMat(rowIdx,:) = HorzDisplacement(:)';
        slopeIncline(rowIdx,:) = side(:)';
        
        %-- update phedi location using the relative shif of the fit
%         NewDisplacement = HorzDisplacementMat(rowIdx,:);
%         PhediCoordinatesUpdate = HorzDisplacementMat(rowIdx,:);
        PhediLocation(rowIdx+1,:) = HorzDisplacementMat(rowIdx,:);
    end
    PhediLocation = PhediLocation(2:end,:);
end
