%[residualsCell, RMS_vec, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStruct, varargin)
%varargin: [reduceSin, plotFits, plotInBlock, XLim, Handle]
% defaults: (1,0,1,[-0.05 -0.01])
%
% phedi_UxBestFitVerticalShift will calculate the shift needed for a good
% fit between LEFM prediction and the Ux measurement.
%


function [residualsCell, RMS_vec, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStruct, varargin)
    
    
    %% set defaults
    [reduceSin, plotFits, plotInBlock, XLim, Handle] = setDefaults4function(varargin,1,0,1,[-0.05 -0.01],[]);
    
    %% set variables
    [Cd, Cs, Cr]=CrackSolutionMaterialProperties;
    PhediData = PhediStruct.PhediData;
    spaceAxis = PhediData.x_mins_x_tip;
    PhediLocation = PhediData.PhediLocation;
    PhotoLocation = PhediStruct.PhediData.PhotoLocation;
    UpperBlockLength = PhediStruct.ExperimentData.UpperBlockLength*1e-3;
    Nphedis = size(PhediLocation,2);
    Phedis2Zero = bsxfun(@minus,PhediLocation,mean(PhediLocation(1:30,:),'omitnan'));
    if reduceSin
        PhediLocReduced = phedi_reduceSinFromLocation(PhediData);
        PhediLocReduced = bsxfun(@minus,PhediLocReduced,mean(PhediLocReduced(1:30,:),'omitnan'));
    else
        PhediLocReduced = Phedis2Zero;
    end
    
    %-- set limits for relevant plotting
    if plotInBlock && isfield(PhediData,'inBlockLogicals')
        inBlockLogicals = logical(PhediData.inBlockLogicals);
    else
        inBlockLogicals = true(size(spaceAxis));
    end
    
    %% pick plot handle
    if plotFits
        if isempty(Handle)
            figure;
        else
            if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
                H = figure(Handle);
                H.Name = ['Event ',num2str(Event),' RowOver w/ Phedis'];
            elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
                axes(Handle);
            else
                error('''Handle'' isn''t axes nor figure');
            end
        end
        
        hold on;
        Colors = MyVaryColor(size(PhediLocReduced,2));
        %     Colors = distinguishable_colors(size(PhediLocReduced,2));
    end
    
    %% define region for fitting
    %--- find rebound from elastic waves, emitted from crack tip, returned
    %-- from boundary.
    travelDist = 2*(UpperBlockLength-PhotoLocation);
    travelTimeRalgh = travelDist/Cr;
    travelTimeShear = travelDist/Cs;
    travelTimeLongtdl = travelDist/Cd;
    timeAxis = PhediData.t_mins_t_tip;
    [~,I] = min(abs(timeAxis - travelTimeRalgh),[],1);
    Ilinear = sub2ind(size(timeAxis),I,1:Nphedis);
    returnWaveLoc = mean(spaceAxis(Ilinear));
    blockEdgeLoc = PhotoLocation-UpperBlockLength+0.01;
    regionEdgeLoc = max([returnWaveLoc,blockEdgeLoc]);
    if regionEdgeLoc>=max(XLim)
        regionEdgeLoc = max(XLim)-0.01;
    end
    XLim = [regionEdgeLoc, max(XLim)];
    
    %% best fit
    residualsCell = cell(size(PhediLocReduced,2),1);
    RMS_vec = zeros(size(PhediLocReduced,2),1);
    shifts_vec = zeros(size(PhediLocReduced,2),1);
    
    for i = 1:size(PhediLocReduced,2)
        inRegion = spaceAxis(:,i)<=max(XLim) & spaceAxis(:,i)>=min(XLim);
        RegionLogic = inBlockLogicals(:,i) & inRegion;
        XX = spaceAxis(RegionLogic,i);
        YY = PhediLocReduced(RegionLogic,i);
        Y_model = interp1(PhediStruct.solAtInter.x,PhediStruct.solAtInter.Ux,XX,'v5cubic');
        
        %-- shift for first two times
        shiftsVec = -peak2peak(Y_model):1e-7:peak2peak(Y_model);
        ShiftedsMat = bsxfun(@minus,Y_model(:),shiftsVec);
        DiffMat = bsxfun(@minus,ShiftedsMat,YY(:));
        [minRMS,I] = min(rms(DiffMat,1));
        minShift = shiftsVec(I);
        residualsCell{i} = DiffMat(:,I);
        if isempty(minRMS); minRMS=nan; end;
        if isempty(minShift); minShift=nan; end;
        RMS_vec(i) = minRMS;
        shifts_vec(i) = minShift;
        
        if plotFits
            plot(spaceAxis(:,i),PhediLocReduced(:,i),...
                'Marker','s','MarkerSize',3,'MarkerFaceColor',[0.6 0.6 0.6],...
                'LineWidth',1,'Color',Colors(i,:));
            plot(PhediStruct.solAtInter.x,PhediStruct.solAtInter.Ux-minShift,'-','Color',Colors(i,:),'LineWidth',2);
        end
    end
    
    
    
end
