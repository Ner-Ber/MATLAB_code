function name = wavesPlotStatistics(linesCell, PhediStructCell,Color)

%% create incline cells
inclines = {};
for i = 1:length(linesCell)
    if ~isempty(linesCell{i})
        inclines{i} = squeeze(2e-4*diff(linesCell{i}(:,2,:),1,1)./diff(linesCell{i}(:,1,:),1,1));
    end
end

%% organize inclines to matrixes:
inclinesMat = nan(max(cellfun(@length, inclines)),length(linesCell));
for i = 1:length(inclines)
    if ~isempty(inclines{i})
        inclinesMat(1:length(inclines{i}),i) = inclines{i}(:);
    end
end

%% create F Vectors
F_vec = zeros(1,length(linesCell));
for i = 1:length(PhediStructCell)
    try
        F_vec(i) = mean(PhediStructCell{i}.SgData.F);
    catch
    end
end

%% create N Vectors
N_vec = zeros(1,length(linesCell));
for i = 1:length(PhediStructCell)
    try
        N_vec(i) = mean(PhediStructCell{i}.SgData.N);
    catch
    end
end

%% plot boxplots
%--- prepare data to plot
meanVals = mean(abs(inclinesMat),'omitnan');
upperBound = abs(meanVals-max(abs(inclinesMat),[],1,'omitnan'));
lowerBound = abs(meanVals-min(abs(inclinesMat),[],1,'omitnan'));
%--- name
name = ['N=',num2str(mean(N_vec(3:end)))];
%--- plot
figure(gcf);
errorbar(F_vec, meanVals,lowerBound,upperBound,'*','Color',Color);


% figure; 
% boxplot(abs(inclinesMat),'Whisker',inf);


end