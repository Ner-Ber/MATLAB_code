%% set data
DataStruct = PhediStruct_3_11;
EvNum = 1;
preKink = -4.6e-4;
myfitStruct = phedi_fitSineToPhediTraj(DataStruct{EvNum}.PhediData,[],2,1);
%-- parameters:
t_min = 0.005/DataStruct{EvNum}.BigPicRotStruct.AvgCf_scratches;
t_max = 0.03/DataStruct{EvNum}.BigPicRotStruct.AvgCf_scratches;
Shifts = linspace(-3*1e-3,3*1e-3,32);
selectedPhedis = [4 6 8 12 14 16 23];

%% create fixed matrix
PhediLoc = DataStruct{EvNum}.PhediData.PhediLocation;
t_mins_t_tips = bsxfun(@minus,DataStruct{EvNum}.PhediData.timeVec,DataStruct{EvNum}.PhediData.t_tips(:)');

N = selectedPhedis; %size(t_mins_t_tips,2);
PhediLocFixed = nan(size(PhediLoc));
for i=N
    PhediLocFixed(:,i) = PhediLoc(:,i)-myfitStruct.model(myfitStruct.b_param{i},myfitStruct.Shifts(i),t_mins_t_tips(:,i));
end

%% plot shifted figures;
MDLs = {};
RMSE_cell = cell(length(Shifts),1);
mean_RMSE = [];
std_RMSE = [];
% N = size(PhediLocFixed,2);
FigColors = MyVaryColor(size(PhediLocFixed,2));
n = floor(sqrt(length(Shifts)));
m = ceil(length(Shifts)/n);
figure;
ha = tight_subplot(n,m,[.02 .03],[.1 .01],[.01 .01]);
for i = 1:length(Shifts)
    MDLs{i} = {};
    RMSE_vec = nan(length(N),1);
    %     figure;
    axes(ha(i));
    hold on;
    for j = N
        %--- plot
        LOGIC = t_mins_t_tips(:,j)>t_min & t_mins_t_tips(:,j)<t_max;
        hold all;
        loglog(t_mins_t_tips(LOGIC,j)+Shifts(i),PhediLocFixed(LOGIC,j),'.-','Color',FigColors(j,:));
        
        %--- fit linear
        %         [P,S] = polyfit(log10(t_mins_t_tips(LOGIC,j)+Shifts(i)),log10(PhediLocFixed(LOGIC,j)),1);
        MDLs{i}{j} = fitlm(real(log10(t_mins_t_tips(LOGIC,j)+Shifts(i))),...
            real(log10(PhediLocFixed(LOGIC,j))));
        RMSE_vec(j) = MDLs{i}{j}.Rsquared.Adjusted;
    end
    RMSE_cell{i} = RMSE_vec;
    mean_RMSE(i) = mean(RMSE_vec);
    std_RMSE(i) = std(RMSE_vec);
    hold all;
    x = linspace(t_min,t_max,50)+Shifts(i);
    y = sqrt(x);
    Inter_i = mean(cellfun(@(A) A.Coefficients{1,1},MDLs{i}(~cellfun(@isempty,MDLs{i}))));
    x1_i = mean(cellfun(@(A) A.Coefficients{2,1},MDLs{i}(~cellfun(@isempty,MDLs{i}))));
    loglog(x,.001*y,'k-o','LineWidth',1.5);
    loglog(x,(10^Inter_i)*(x).^x1_i,'c-o','LineWidth',1.5);
    title(['shift=',num2str(Shifts(i))]);
    set(gca,'XScale','log','YScale','log')
end

%% plot fits statistics
figure; plot(Shifts,mean_RMSE,'*');
xlabel('shift'); ylabel('mean RMSE');
figure; plot(Shifts,std_RMSE,'*');
xlabel('shift'); ylabel('std RMSE');

Close2Half = cellfun(@(B) mean(cellfun(@(A) A.Coefficients{2,1},B))-0.5,MDLs);
figure; plot(Shifts,Close2Half,'*');
xlabel('shift'); ylabel('power-0.5');

%% fit parabola and find optimal shift
% [~,M] = min(mean_RMSE);
% p = polyfit(Shifts(M-2:M+2),mean_RMSE(M-2:M+2),2);
% minShift = -p(2)/(2*p(1));
minShift = find_x0_at_y0(0,Close2Half,Shifts);

%% build model for optimal shift
intercept = [];
x1 = [];
mdl = {};
for j = N
    %--- plot
    LOGIC = t_mins_t_tips(:,j)>t_min & t_mins_t_tips(:,j)<t_max;
    %--- fit linear
    mdl{j} = fitlm(log10(t_mins_t_tips(LOGIC,j)+minShift),log10(PhediLocFixed(LOGIC,j)));
    intercept(j) = mdl{j}.Coefficients{1,1};
    x1(j) = mdl{j}.Coefficients{2,1};
end
interceptMean = mean(intercept);
x1Mean = mean(abs(x1));

%% plot fixed data wit sqrt
figure;
xx = linspace(0,max(t_mins_t_tips(:)),500);
yy_sqrt = [zeros(1,499),(10^interceptMean)*sqrt(xx)];
yy_fit = [zeros(1,499),(10^interceptMean)*((xx).^x1Mean)];
xx = [-fliplr(xx(2:end)),xx];
% N = size(t_mins_t_tips,2);
% FigColors = MyVaryColor(N);
Markers = repmat({'*'},1,N);
Markers(DataStruct{EvNum}.PhediData.slopeIncline==-1) = {'.'};
Markers(DataStruct{EvNum}.PhediData.slopeIncline==0) = {'diamond'};

for i=N
    hold on;
    plot(t_mins_t_tips(:,i),PhediLocFixed(:,i),...
        '.-','Color',FigColors(i,:),'Marker',Markers{i});
end
hold on;
plot(xx-minShift,yy_fit,'k','LineWidth',1.3);
plot(xx-minShift,yy_sqrt,'Color',rgb('RosyBrown'),'LineWidth',1.3);