exper = my_dir
moveOn = 1;
ExperNum = 1;
EvNum=1;
failNum=0;
BigPicRowOverTimeStruct={};
while moveOn
    try
%         BigPicRowOverTimeStruct{EvNum} = IDT_PlotRowOverTime(exper{1},EvNum,[],[],[],[],[],[],[],0);
%         BigPicRowOverTimeStruct{EvNum} = IDT_fixRotStructXaxis(BigPicRowOverTimeStruct{EvNum});
%         figure;
%         IDT_PlotRowOverTime(BigPicRowOverTimeStruct{EvNum});
        
        figure;
        subplot(2,1,1);
%         BigPicRowOverTimeStruct{EvNum} = IDT_PlotRowOverTime(exper{ExperNum},EvNum);
        BigPicRowOverTimeStruct{EvNum} = phantom_BigPicStruct(exper{ExperNum},EvNum,3e-3,1e-10,3e-3,1,1,'all');
        IDT_PlotRowOverTime(BigPicRowOverTimeStruct{EvNum});
        
        scratchRegionMeters = [0.07 0.09];
        xq = linspace(scratchRegionMeters(1),scratchRegionMeters(2),50);
        scartchLocationsInPix = interp1(BigPicRowOverTimeStruct{EvNum}.x,1:length(BigPicRowOverTimeStruct{EvNum}.x),xq);
        scartchLocationsInPix = unique([scartchLocationsInPix,round(min(scartchLocationsInPix)):round(max(scartchLocationsInPix))]);
        frontTimeInterped = interp1(BigPicRowOverTimeStruct{EvNum}.frontStepsPix,BigPicRowOverTimeStruct{EvNum}.frontStepsTime,scartchLocationsInPix);
        MeanFronVel = mean(diff(BigPicRowOverTimeStruct{EvNum}.x))*mean(diff(scartchLocationsInPix)./diff(frontTimeInterped))./mean(diff(BigPicRowOverTimeStruct{EvNum}.t));
        h1=get(gca,'title');
        title({h1.String,['MeanFronVel= ',num2str(MeanFronVel)]});
        
        subplot(2,1,2);
        phantom_PlotRowOverTime(exper{ExperNum},EvNum);
        
        EvNum = EvNum+1;
     catch
        EvNum = EvNum+1;
        failNum=failNum+1;
        if failNum>=5
            break
        end
    end
    
end