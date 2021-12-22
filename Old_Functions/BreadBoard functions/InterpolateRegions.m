function [currentXgridI,currentYgridI,interpped_frames,Extremas] = ...
    InterpolateRegions(Xgrid,Ygrid,InterpFactor,InterpMethod,frame_range,...
    coordinatesCell,fixed_frames_cr,regionIdx,Extremas)
    
x_1 = min(coordinatesCell{regionIdx}([1,3]));
x_2 = max(coordinatesCell{regionIdx}([1,3]));
y_1 = min(coordinatesCell{regionIdx}([2,4]));
y_2 = max(coordinatesCell{regionIdx}([2,4]));
currentRegionFrames = fixed_frames_cr(y_1:y_2,x_1:x_2,:);

%---interpolate frames
    currentXgrid = Xgrid(y_1:y_2,x_1:x_2);
    currentYgrid = Ygrid(y_1:y_2,x_1:x_2);
    [currentXgridI, currentYgridI] = meshgrid(linspace(x_1,x_2,(x_2-x_1)*InterpFactor),...
    linspace(y_1,y_2,(y_2-y_1)*InterpFactor));

interpped_frames = zeros(size(currentXgridI,1),size(currentXgridI,2),length(frame_range));
    for timeStep = 1:length(frame_range)
        interpped_frames(:,:,timeStep) = interp2(currentXgrid, currentYgrid,...
            currentRegionFrames(:,:,timeStep), currentXgridI, currentYgridI, InterpMethod);
        %---find extremas
        [frameMAX,IMAX,~,~] = extrema2(interpped_frames(:,:,timeStep));
        Extremas{regionIdx,timeStep} = [frameMAX,IMAX];
    end
end