
%% Makes matrixes of the mass center at each frame and each region
Xcenter = zeros(length(coordinatesCell),length(frame_range));
Ycenter = zeros(length(coordinatesCell),length(frame_range));
 for regionIdx = 1:length(coordinatesCell)
     for timeStep = 1:length(frame_range)
         [X_c,Y_c]=Center_of_Mass(interppesFramesCell{regionIdx}(:,:,timeStep),InterpFactor);
         Xcenter(regionIdx,timeStep)=X_c+coordinatesCell{regionIdx}(1)-1;
         Ycenter(regionIdx,timeStep)=Y_c+coordinatesCell{regionIdx}(2)-1;
     end
 end
 
%% plot Xcenter
figure
hold all
for regionIdx=1:length(coordinatesCell)
    plot(frame_range,Xcenter(regionIdx,:),'.')
end
hold off

%% plot Ycenter
figure
hold all
for regionIdx=1:length(coordinatesCell)
    plot(frame_range,Ycenter(regionIdx,:),'.')
end
hold off