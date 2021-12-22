function [vallyPairsCell,asper_num] = Movie_define_asperities_for_scratches(RowOverTimeCell)

vallyPairsCell = cell(size(RowOverTimeCell));
asper_num = zeros(1,length(RowOverTimeCell));
for i = 1:length(RowOverTimeCell)
    signal = RowOverTimeCell{i}(1,:);
    
    N = 80;
    smoothFactor = round(length(signal)/N);
    SmoothedSignal = smooth(signal,smoothFactor);
    SmoothedSignal = SmoothedSignal-(mean(SmoothedSignal)*0.5);
    
    trenchesLogical = SmoothedSignal>0;
    lable = bwlabel(trenchesLogical);
    trenchesLogical(lable==1|lable==max(lable)) = 0;
    Dtrenches = [0;diff(trenchesLogical)];
    LeftMark = find(Dtrenches==1);
    RightMark = find(Dtrenches==-1);
    vallyPairsCell{i} = [LeftMark(:),RightMark(:)];
    asper_num(i) = length(LeftMark);
end
end