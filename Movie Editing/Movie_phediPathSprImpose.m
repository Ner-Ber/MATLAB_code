function Movie_phediPathSprImpose(analyzePhediCell, analyzeAsperityCell, Event ,phediNum)
%Movie_phediPathSprImpose(analyzePhediCell, analyzeAsperityCell, Event ,phediNum)
% thie function will plot the super position of signals (from RowOverTime)
% in a form of a certain specified phedi anchored to the same location.

relevantFrames = analyzePhediCell{Event}.firstFrame4Phedi:analyzePhediCell{Event}.lastFrame4Phedi;
frameIdx = 1:length(relevantFrames);
Colors = MyVaryColor(length(relevantFrames));

basicX = analyzeAsperityCell{Event}.spatialVec;

figure;
hold on;
for i = frameIdx
    Strc_Fctr = analyzePhediCell{Event}.StrechingFactorsMat(i,phediNum);
    thisX = basicX-(analyzePhediCell{Event}.measuredPhedisFromPlot(phediNum)-analyzePhediCell{Event}.PhediLocation(i,phediNum));
    thisY = analyzeAsperityCell{Event}.RowOverTime(relevantFrames(i),:);
    plot(thisX, Strc_Fctr*thisY,'Color',Colors(i,:))
    title(['Event=',num2str(Event),'  phediNum=',num2str(phediNum),...
        '  initial loc=',num2str(analyzePhediCell{Event}.measuredPhedisFromPlot(phediNum))]);
end


end