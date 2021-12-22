function [vallyPairsCell,asper_num] = Movie_define_asperities_for_scratches(RowOverTime)


signal = RowOverTime(1,:);

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
vallyPairsCell = [LeftMark(:),RightMark(:)];
asper_num = length(LeftMark);

end