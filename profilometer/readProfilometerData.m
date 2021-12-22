
mainFolder = 'G:\Frics\2016-3-22\';
clear fileName
folderName = '\20160322 steel smooth sample to taste the stage\';
d = dir([mainFolder , folderName, '\']);
d(1:2) = [];
filesNameAll = {d.name};
result = cellfun(@(s) length(s)==23 , filesNameAll);
fileNames = filesNameAll(result);


for indFileName = 1:length(fileNames)
fileNameTemp = fileNames{indFileName};
    
mainData = importdata([mainFolder , folderName , fileNameTemp]);    
data = mainData.data;

clear middleInd tempNewX
% data(:,1) is encoder 1, not sure what it is
% data(:,2) is encoder 2 which is the distance from the block (z
% location) - the depth
% data(:,3) is x location - the dirction parallel to the table
% data(:,4) is y location - the dirction perpendicular to the table
% data(:,5) is the nubmer of measurement (of row).

%% Plot Altitude Vs location
% USing the approximate inclination of the line to reduce all the points which are the edges or errors,
% by assuming if it is too far from the inclined line it is not par of the
% data

figure;
hold all;
for i = 1:max(data(:,5))  
    plot(data(data(:,5)==i,2));
end
title(sprintf(' %s distnace measurement of sample' , [fileNameTemp(1:10) , ' ' ,fileNameTemp(12:19)]));
xlabel('measurement points');
ylabel('altitude (\mum)');
setFigFontSize(gcf,16,18);

figure;
for indRow = 1:max(data(:,5))
    hold all;
    % find the limits of the measurement area:
    % We ignore from all the points which are far eonough (around 50 micron) from the mean inclination
    % line, and then do again linear fit by detrend
    depth = data(data(:,5)==indRow,2);
    X = data(data(:,5)==indRow,3);
    newDepth = depth;
    newX = X;
    badPointsZero = newDepth==0;
    newDepth(badPointsZero) = [];
    newX(badPointsZero) = [];
    middleInd = length(newDepth)/2; %taking only the middle part to approximately evaluate the slope of the surface%
    tempNewDepth = newDepth(max(1,middleInd - 1000):min(length(newDepth),middleInd + 1000));
    tempNewX = newX(max(1,middleInd - 1000):min(length(newDepth),middleInd + 1000));
    p = polyfit(tempNewX,tempNewDepth,1); %polyfit to linear
    testOutsideLinear = 100; %in micron. The test is if a point is distant from the line more than testOutsideLinear.
    badPoints = abs(newDepth - (p(1)*newX+p(2))) > testOutsideLinear;
    newDepth(badPoints) = [];
    newX(badPoints) = [];
    y = detrend(newDepth);
    hold all;
    plot(newX/1000 , y , '.-'); %/1000 to be in mm

end
title(sprintf(' %s roughness of sample' , [fileNameTemp(1:10) , ' ' ,fileNameTemp(12:19)]));
xlabel('X (mm)');
ylabel('altitude (\mum)');
setFigFontSize(gcf,16,18);


end