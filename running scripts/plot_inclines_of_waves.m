%% create incline cells
%--- exp 3
inclines_1to7_Exp3 = {};
for i = 1:length(linesCell_1to7_Exp3)
    if ~isempty(linesCell_1to7_Exp3{i})
        inclines_1to7_Exp3{i} = squeeze(2e-4*diff(linesCell_1to7_Exp3{i}(:,2,:),1,1)./diff(linesCell_1to7_Exp3{i}(:,1,:),1,1));
    end
end
%--- exp 4
inclines_1to30_Exp4= {};
for i = 1:length(linesCell_Exp4)
    if ~isempty(linesCell_Exp4{i})
        inclines_1to30_Exp4{i} = squeeze(2e-4*diff(linesCell_Exp4{i}(:,2,:),1,1)./diff(linesCell_Exp4{i}(:,1,:),1,1));
    end
end
%--- exp 4
inclines_1to28_Exp5= {};
for i = 1:length(linesCell_Exp5)
    if ~isempty(linesCell_Exp5{i})
        inclines_1to28_Exp5{i} = squeeze(2e-4*diff(linesCell_Exp5{i}(:,2,:),1,1)./diff(linesCell_Exp5{i}(:,1,:),1,1));
    end
end

%% organize inclines to matrixes:
%--- exp 3
inclines_Exp3_mat = nan(max(cellfun(@length, inclines_1to7_Exp3)),length(linesCell_1to7_Exp3));
for i = 1:length(inclines_1to7_Exp3)
    if ~isempty(inclines_1to7_Exp3{i})
        inclines_Exp3_mat(1:length(inclines_1to7_Exp3{i}),i) = inclines_1to7_Exp3{i}(:);
    end
end

%--- exp 4
inclines_Exp4_mat = nan(max(cellfun(@length, inclines_1to30_Exp4)),length(linesCell_Exp4));
for i = 1:length(inclines_1to30_Exp4)
    if ~isempty(inclines_1to30_Exp4{i})
        inclines_Exp4_mat(1:length(inclines_1to30_Exp4{i}),i) = inclines_1to30_Exp4{i}(:);
    end
end

%--- exp 5
inclines_Exp5_mat = nan(max(cellfun(@length, inclines_1to28_Exp5)),length(linesCell_Exp5));
for i = 1:length(inclines_1to28_Exp5)
    if ~isempty(inclines_1to28_Exp5{i})
        inclines_Exp5_mat(1:length(inclines_1to28_Exp5{i}),i) = inclines_1to28_Exp5{i}(:);
    end
end

%% create F Vectors
F_vec_Exp3 = zeros(1,length(linesCell_1to7_Exp3));
for i = 1:length(PhediStructCell_1to7_Exp3)
    try
        F_vec_Exp3(i) = mean(PhediStructCell_1to7_Exp3{i}.SgData.F);
    catch
    end
end

F_vec_Exp4 = zeros(1,length(linesCell_Exp4));
for i = 1:length(PhediStructures_1to30_Exp4)
    try
        F_vec_Exp4(i) = mean(PhediStructures_1to30_Exp4{i}.SgData.F);
    catch
    end
end

F_vec_Exp5 = zeros(1,length(linesCell_Exp5));
for i = 1:length(PhediStructures_1to28_Exp5)
    try
        F_vec_Exp5(i) = mean(PhediStructures_1to28_Exp5{i}.SgData.F);
    catch
    end
end

%% create N Vectors
N_vec_Exp3 = zeros(1,length(linesCell_1to7_Exp3));
for i = 1:length(PhediStructCell_1to7_Exp3)
    try
        N_vec_Exp3(i) = mean(PhediStructCell_1to7_Exp3{i}.SgData.N);
    catch
    end
end

N_vec_Exp4 = zeros(1,length(linesCell_Exp4));
for i = 1:length(PhediStructures_1to30_Exp4)
    try
        N_vec_Exp4(i) = mean(PhediStructures_1to30_Exp4{i}.SgData.N);
    catch
    end
end

N_vec_Exp5 = zeros(1,length(linesCell_Exp5));
for i = 1:length(PhediStructures_1to28_Exp5)
    try
        N_vec_Exp5(i) = mean(PhediStructures_1to28_Exp5{i}.SgData.N);
    catch
    end
end



%% plot boxplots

figure;
hold on;
plot(F_vec_Exp3, mean(abs(inclines_Exp3_mat),'omitnan'),'*');
plot(F_vec_Exp4, mean(abs(inclines_Exp4_mat),'omitnan'),'*');
plot(F_vec_Exp5, mean(abs(inclines_Exp5_mat),'omitnan'),'*');


figure; 
boxplot(abs(inclines_Exp3_mat),'Whisker',inf);

figure; 
boxplot(abs(inclines_Exp4_mat),'Whisker',inf);

figure; 
boxplot(abs(inclines_Exp5_mat),'Whisker',inf);