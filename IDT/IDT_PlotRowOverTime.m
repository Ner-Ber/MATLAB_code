function RowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,varargin)
%RowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,pre_time,intervalTime,post_time,smt,lineNum, frontMarkLocation, figHandle, ifPlot)
%will plot the row over time of the IDT camera
%
%exp_dir - path of the experiment folder
%eventNum - event number
%VARARGIN (DEFAULTS):
%pre_time - time before event to display (def = 3e-3 (s))
%intervalTime - interval in display frames (def = 1e-6 (s))
%post_time - time after event to display (def = 3e-3 (s))
%smt - smooth factor along the x axis (default=0)
%lineNum - lines to consider (default = 'all');
%frontMarkLocation - region in which to mark front upon figure (in meters) (default=[0 0.15])
%figHandle - handle to plot figure (default=gcf)
%ifPlot -  option to plot (default=1, yes)
%
%   !!!!!!!::
% if instead of exp_dir a structure identical in fields as the output is
% inserted, the function will plot the data from the structure.
% In this case:
%   -) if there is no second input, IDT_PlotRowOverTime will plot the
% regular heat-map
%   -) if the seconf input is 'I' IDT_PlotRowOverTime will plot the
%   not-normalized image.
%   -) if the second input is a positive integer IDT_PlotRowOverTime will
%   normalize by that row
%   -) if the second input is in the bounds of the time vector
%   IDT_PlotRowOverTime will normalize bt that specific time.


if isstruct(exp_dir)
    %     eventNum = nan;
    plotFromStruct = 1;
else
    plotFromStruct = 0;
end


G = groot;
if isempty(G.CurrentFigure)
    [pre_time,intervalTime,post_time,smt,lineNum, frontMarkLocation, figHandle, ifPlot] =...
        setDefaults4function(varargin,...
        3e-3,1e-6,3e-3,0,'all',[0 0.15], [], 1);
else
    [pre_time,intervalTime,post_time,smt,lineNum, frontMarkLocation, figHandle, ifPlot] =...
        setDefaults4function(varargin,...
        3e-3,1e-6,3e-3,0,'all', [0 0.15], gcf, 1);
end

if ~plotFromStruct
    
    ims = IDT_imagesWithData(exp_dir,eventNum,pre_time,intervalTime,post_time,smt,lineNum);
    
    firRowMat = repmat(ims.lines(1,:),size(ims.lines,1),1);
    NormalizedRowOverTime = ims.lines./firRowMat;
    
%     frontStructure = IDT_FindCfFromRowOverTime(ims.lines);
    frontStructure = IDT_FindCfFromRowOverTimeControlFallingPercnt(ims.lines,[],2);
    frontIntepStruct = IDT_interpolateFrontBetweenSteps(frontStructure);
    
    %--- get velocity related vecs in [m/s] and [m]
    frontVelLoc_interpM = interp1(1:length(frontStructure.fronLoc),ims.x,frontIntepStruct.VelLoc,'linear','extrap');
    frontVel_interpMperS = frontIntepStruct.FrontVel*ims.fps*ims.res;
    
    %--- find average velocity in scratch region
    scratchRegionMeters = [0.07 0.09];
    xq = linspace(scratchRegionMeters(1),scratchRegionMeters(2),50);
    scartchLocationsInPix = interp1(ims.x,1:length(ims.x),xq);
    scartchLocationsInPix = unique([scartchLocationsInPix,round(min(scartchLocationsInPix)):round(max(scartchLocationsInPix))]);
    frontTimeInterped = interp1(frontStructure.StepsPix,frontStructure.StepTime,scartchLocationsInPix);
    AvgCf_scratches = mean(diff(ims.x))*mean(diff(scartchLocationsInPix)./diff(frontTimeInterped))./mean(diff(ims.t));
    
    if nargout>0
        
        RowOverTimeStruct = struct;
        RowOverTimeStruct.DataMat = ims.lines;
        RowOverTimeStruct.DataMatNorm = NormalizedRowOverTime;
        RowOverTimeStruct.res = ims.res;
        RowOverTimeStruct.fps = ims.fps;
        RowOverTimeStruct.x = ims.x;
        RowOverTimeStruct.t = ims.t;
        RowOverTimeStruct.details = [ims.Date,'   ',ims.Exp,'   event=',num2str(ims.EventNum)];
        RowOverTimeStruct.xlabel = 'x [m]';
        RowOverTimeStruct.ylabel = 't [s]';
        RowOverTimeStruct.frontLoc = frontStructure.fronLoc;
        RowOverTimeStruct.frontStepsPix = frontStructure.StepsPix;
        RowOverTimeStruct.frontStepsTime = frontStructure.StepTime;
        RowOverTimeStruct.frontVel = frontStructure.FrontVel;
        RowOverTimeStruct.frontVelLoc = frontStructure.VelLoc;
        RowOverTimeStruct.frontLoc_interp = frontIntepStruct.xq;
        %     RowOverTimeStruct.frontTime_interp = frontIntepStruct.tq;
        RowOverTimeStruct.frontTime_interp = frontIntepStruct.tq+ims.t(1)*ims.fps-1;    % time vector for front in units of frames
        RowOverTimeStruct.frontTimeOriginal_interp = frontIntepStruct.tq_original+ims.t(1)*ims.fps-1;    % same as frontTime_interp but not smoothed
        RowOverTimeStruct.frontVel_interp = frontIntepStruct.FrontVel;
        RowOverTimeStruct.frontVelLoc_interp = frontIntepStruct.VelLoc;
        RowOverTimeStruct.frontVelLoc_interpM = frontVelLoc_interpM;
        RowOverTimeStruct.frontVel_interpMperS = frontVel_interpMperS;
        RowOverTimeStruct.AvgCf_scratches = AvgCf_scratches;
        RowOverTimeStruct.BlockEdgePixles = ims.BlockEdgePixles;
    end
    
    if ifPlot
        if ~isempty(findobj('type','figure')) && isa(figHandle,'matlab.graphics.axis.Axes')
            axes(figHandle)
        elseif isempty(findobj('type','figure'));
            figure;
        else
            figure(figHandle);
        end
        imagesc([ims.x(1) ims.x(end)],[ims.t(1) ims.t(end)],NormalizedRowOverTime);
        title([ims.Date,'   ',ims.Exp,'   event=',num2str(ims.EventNum)]);
        xlabel('x [m]'); ylabel('t [s]');
        colorbar;
        ax = gca;
        ax.YDir = 'normal';
        caxis([0.5 1.1]);
        
        hold on;
        %--- plot front interpolated
        markRegionLogical = ims.x>=min(frontMarkLocation)& ims.x<=max(frontMarkLocation);
        x_plot = ims.res*frontIntepStruct.xq+ims.x(1);
        t_plot = (frontIntepStruct.tq+ims.t(1)*ims.fps-1)/ims.fps;
        plot(x_plot(markRegionLogical), t_plot(markRegionLogical), 'r','LineWidth',1.5);
        
        %--- plot front steps
        x_plot = ims.x(frontStructure.StepsPix);
        x_plot(x_plot<min(frontMarkLocation)| x_plot>max(frontMarkLocation)) = nan;
        t_plot = ims.t(frontStructure.StepTime);
        plot(x_plot, t_plot, 'b*','LineWidth',1.5);
        
        hold off;
    end
else
    RotStruct = exp_dir;
    plotIntensity = 0;
    %-- plot normalized or intensity?
    if exist('eventNum','var') && ~isempty(eventNum)
        if ischar(eventNum)
            plotIntensity = (strcmpi(eventNum,'Intensity') || strcmpi(eventNum,'Intens') || strcmpi(eventNum,'I'));
            if ~plotIntensity
                relevantROT = RotStruct.DataMatNorm;
            else
                relevantROT = RotStruct.DataMat;
            end
        elseif isnumeric(eventNum) %-- if isnumeris the DataMat will be normalized by time or row number
            if eventNum>=max(RotStruct.t) && mod(eventNum,1)==0 %-- this is the case where 'eventNum' represents row number
                normRow = eventNum;
            else    %-- this is the case where 'eventNum' represents row a time
                [~,normRow] = min(abs(RotStruct.t-eventNum));
            end
            relevantROT = bsxfun(@rdivide,RotStruct.DataMat,RotStruct.DataMat(normRow,:));
        end
    elseif ~exist('eventNum') || isempty(eventNum)
        relevantROT = RotStruct.DataMatNorm;
    end
    
    if ~isempty(findobj('type','figure')) && isa(figHandle,'matlab.graphics.axis.Axes')
        axes(figHandle)
    elseif isempty(findobj('type','figure'));
        figure;
    else
        figure(figHandle);
    end
%     
%     if isempty(findobj('type','figure'));
%         figure;
%     else
%         figure(figHandle);
%     end
    imagesc([RotStruct.x(1) RotStruct.x(end)],[RotStruct.t(1) RotStruct.t(end)],relevantROT);
    title({RotStruct.details,['mean Cf @ pattern=',num2str(round(RotStruct.AvgCf_scratches))]});
    xlabel('x [m]'); ylabel('t [s]');
    colorbar;
    ax = gca;
    ax.YDir = 'normal';
    if ~plotIntensity
        caxis([0.5 1.1]);
    end
    
    hold on;
    %--- plot intepolated front
    markRegionLogical = RotStruct.x>=min(frontMarkLocation)& RotStruct.x<=max(frontMarkLocation);
    %     x_plot = RotStruct.res*RotStruct.frontLoc_interp+RotStruct.x(1);
    x_plot = RotStruct.x;
    t_plot = (RotStruct.frontTime_interp)/RotStruct.fps;
    plot(x_plot(markRegionLogical), t_plot(markRegionLogical), 'r','LineWidth',1.5);
    %--- plot front steps
    x_plot = RotStruct.x(RotStruct.frontStepsPix);
    x_plot(x_plot<min(frontMarkLocation)| x_plot>max(frontMarkLocation)) = nan;
    t_plot = RotStruct.t(RotStruct.frontStepsTime);
    plot(x_plot, t_plot, 'b*','LineWidth',1.5);
    
    hold off;
end

end