function [LocDeTrended, usdTime, meanVelMeas, meanVelTrndd] = phedi_compMeanField2Measure(timeVecFrms,timeIntrvlFrms,PhediLocPix)
    
    timeLogical = timeVecFrms>=min(timeIntrvlFrms) & timeVecFrms<=max(timeIntrvlFrms);
    
    usdTime = timeVecFrms(timeLogical);
    usdTime = usdTime(:);
    meanVelTrndd = [];
    LocDeTrended = zeros(size(PhediLocPix(timeLogical,:)));
    for i=1:size(PhediLocPix,2)
        p = polyfit(usdTime(:),PhediLocPix(timeLogical,i),1);
        LocDeTrended = PhediLocPix(timeLogical,i)-polyval(p,usdTime);
        meanVelTrndd = cat(1,meanVelTrndd,p(1));
    end
    
    meanVelMeas = mean(bsxfun(@rdivide,diff(PhediLocPix(timeLogical,:),1,1),diff(usdTime)),1);
    
    
end
