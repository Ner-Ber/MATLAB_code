function [multiLineData, SurfPlot, slicesFig] = multipLines_createImagesOfFront2D(multipLinesStruct)
    if ischar(multipLinesStruct)
        multipLinesStruct = load(multipLinesStruct);
        multipLinesStruct = multipLinesStruct.BigPicMultipLines;
    end
    
    %-- save data from all lines
    names = fieldnames(multipLinesStruct);
    multipLinesROT = [];
    x_fronts = [];
    t_fronts = [];
    for n = names(:)'
        thisLineStruct = multipLinesStruct.(n{1});
        %--- save the ROT matrix
        multipLinesROT = cat(3,multipLinesROT,thisLineStruct.DataMatNorm);
        %--- save path trajectory
        x_fronts = cat(1,x_fronts,thisLineStruct.x(:)');
        t_fronts = cat(1,t_fronts, (thisLineStruct.frontTime_interp(:)')/thisLineStruct.fps);
    end
    ax_time = thisLineStruct.t(:)';
    ax_space = thisLineStruct.x(:)';
    
    multiLineData = struct;
    multiLineData.x_fronts = x_fronts;
    multiLineData.t_fronts = t_fronts;
    multiLineData.multipLinesROT = multipLinesROT;
    multiLineData.ax_time = ax_time;
    multiLineData.ax_space = ax_space;
    
    
    %--- plot the ROT on specific places
    x_front_mean = median(x_fronts,1);
    t_front_mean = median(t_fronts,1);
    groove_region = [0.07, 0.09];
    M = abs(bsxfun(@minus, groove_region,x_front_mean(:)));
    [~,groove_region_idx] = min(M,[],1);
    groove_time_bounds = t_front_mean(groove_region_idx);
    N = 7;
    time_samps = linspace(min(groove_time_bounds), max(groove_time_bounds), N);
    
    slicesFig = figure; hold on;
    for i = 1:N
        AX = subplot(N,1,i);
        hold on;
        
        [~,minIdx] = min(abs(ax_time-time_samps(i)));
        snapshot = permute(multipLinesROT(minIdx,:,:),[3,2,1]);
        snap_img = imagesc([ax_space(1), ax_space(end)],[1,size(snapshot,1)], snapshot);
        AX.XTick = [0,0.05,linspace(0.07,0.09,11), 0.1, 0.15];
        if i~=N
            AX.XTickLabel=[];
        end
        xlim([0,0.15]);
        yl = ylabel(sprintf('t=%1.2e [s]', ax_time(minIdx)));
        yl.Rotation=0;
        yl.HorizontalAlignment='right';
        caxis([0.7, 1.1]);
        
        plot([0.07, 0.07],[0.5,size(snapshot,1)+0.5],'r','LineWidth',1);
        plot([0.09, 0.09],[0.5,size(snapshot,1)+0.5],'r','LineWidth',1);
        
        ylim([0.5,size(snapshot,1)+0.5]);
    end
    
    
    %--- plot the 3d surface
    figure;
    hold on;
    L_fronts = bsxfun(@times,ones(1,1280), (1:8)');
%     copper_cmap = colormap(copper);
    SurfPlot = surf(x_fronts, L_fronts, t_fronts);
    SurfPlot.EdgeAlpha=0;
    colorbar
    xlabel('X [m]')
    ylabel('line')
    zlabel('times [s]')
    for i=1:8
        plot3(x_fronts(i,:),i*ones(1,1280), t_fronts(i,:));
    end
    zlim(groove_time_bounds);
    caxis(groove_time_bounds);
    xlim([0.06, 0.1]);
    view(0,0)
    
    
end


