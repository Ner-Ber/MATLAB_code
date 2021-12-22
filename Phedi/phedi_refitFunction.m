function phedi_refitFunction(RowOverTime_normalzd,phediNum,PhediLocationPixPreSmth,firstFrame4Phedi,lastFrame4Phedi)
    
    
    %% parameters
    roundMagntd = 4;
    %%
    [XShiftedMat, RowOverTimeRelevant, phediIntialLocPix] = phedi_reshiftByPhedi(RowOverTime_normalzd,phediNum,PhediLocationPixPreSmth,firstFrame4Phedi,lastFrame4Phedi);
%     SlopeRegions = PhediStruct_7_13.PhediData.SlopeRegions;
    
    
    %% find region to fit tanh to
    relevantRegion = [-4 4];
    relevantLogic = XShiftedMat>=min(relevantRegion) & XShiftedMat<=max(relevantRegion);
    XRelevant = XShiftedMat(relevantLogic);
    Yrelevant = RowOverTimeRelevant(relevantLogic);
    XRelevantRound = round(XRelevant,roundMagntd);
    [Xunique,~,ic] = unique(XRelevantRound);
    Ncount = histc(XRelevantRound, Xunique);
    relevantXcoor = round(min(Xunique):10^-roundMagntd:max(Xunique),roundMagntd);
    NanMat = nan(max(Ncount),length(relevantXcoor));
    for Idx = 1:length(Xunique)
        [~,colIdx] = min(abs(relevantXcoor-mean(XRelevantRound(Idx==ic))));
        NanMat(1:length(Yrelevant(Idx==ic)),colIdx) = Yrelevant(Idx==ic);
        
    end
    %-- choose window for moving variance
    n = max(diff(sort(XRelevantRound)))./10^-roundMagntd;
    N = round(n*1.1);
    N = N + double(~mod(N,2));  % make N odd
    M = max(Ncount);
    kernel = ones(M,N)./(M*N);
    
    A1 = NanMat;
    A1(isnan(NanMat)) = 0;
    A2 = double(~~A1);
    FirstConv = conv2(A1,kernel,'valid');
    ScndConv = conv2(A2,kernel,'valid');
    MovingAvg = FirstConv./ScndConv;
    
    FirstConv = conv2(A1.^2,kernel,'valid');
    ScndConv = conv2(A2.^2,kernel,'valid');
    MovingAvgSqrd = FirstConv./ScndConv;
    
    movingVar = MovingAvgSqrd-MovingAvg.^2;
    X4conv = relevantXcoor(ceil(N/2):(end-floor(N/2)));
    %-- measure variance moving variance
    
    
    
    
end