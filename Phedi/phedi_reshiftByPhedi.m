function [XShiftedMat, RowOverTimeRelevant, phediIntialLocPix] = phedi_reshiftByPhedi(...
        RowOverTime,phediNum,PhediLocationPixPreSmth,firstFrame4Phedi,lastFrame4Phedi,varargin)
    [SkipFact, DoPlot]= setDefaults4function(varargin,5,0);
    
    
    XPreShift = 1:size(RowOverTime,2);
    RowOverTimeRelevant = RowOverTime(firstFrame4Phedi:lastFrame4Phedi,:);
    
    XShiftedMat = bsxfun(@minus,XPreShift,PhediLocationPixPreSmth(:,phediNum));
    phediIntialLocPix = PhediLocationPixPreSmth(1,phediNum);
    
    if DoPlot
        Colors = MyVaryColor(size(RowOverTimeRelevant,1));
        figure; hold on;
        for i=1:SkipFact:size(RowOverTimeRelevant,1)
            plot(XShiftedMat(i,:),RowOverTimeRelevant(i,:),'.','Color',Colors(i,:));
        end
    end
end
