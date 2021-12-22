notEmpty = find(cellfun(@(A) ~isempty(A), PhediStruct));
C_time = cell(size(notEmpty));
Cf = zeros(size(notEmpty));
EvColors = MyVaryColor(length(notEmpty));
figure; hold on;
for i = 1:length(notEmpty)
    try
    EvNum = notEmpty(i);
    [X_kink,Y_kink] = phedi_findCohesiveZoneKink(PhediStruct{EvNum}.PhediData);
    [X_kink2,Y_kink2] = phedi_findCohesiveZoneKink(PhediStruct{EvNum}.PhediData,[1e-5 max(PhediStruct{EvNum}.PhediData.timeVec)-5e-6]);
    C_time{i} = X_kink2-X_kink;
    Cf(i) = PhediStruct{EvNum}.BigPicRotStruct.AvgCf_scratches;
    plot(repmat(Cf(i),1,length(C_time{i})),C_time{i},'o','Color',EvColors(i,:));
    catch
    end
end

xlabel('Cf [m/s]'); ylabel('cohesive time (s)');