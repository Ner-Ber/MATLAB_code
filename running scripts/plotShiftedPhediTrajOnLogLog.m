%% set data
DataStruct = PhediStruct_3_11;
EvNum = 1;
preKink = -4.6e-4;
myfitStruct = phedi_fitSineToPhediTraj(DataStruct{EvNum}.PhediData,[],2,1);
%-- parameters:
t_min = 0.003/DataStruct{EvNum}.BigPicRotStruct.AvgCf_scratches;
t_max = 0.025/DataStruct{EvNum}.BigPicRotStruct.AvgCf_scratches;
Shifts = linspace(-2*1e-5,2*1e-5,20);
selectedPhedis = [4 6 8 12 14 16 23];

%% create fixed matrix
PhediLoc = DataStruct{EvNum}.PhediData.PhediLocation;
t_mins_t_tips = bsxfun(@minus,DataStruct{EvNum}.PhediData.timeVec,DataStruct{EvNum}.PhediData.t_tips(:)');

N = size(t_mins_t_tips,2);
PhediLocFixed = nan(size(PhediLoc));
for i=1:N
    PhediLocFixed(:,i) = PhediLoc(:,i)-myfitStruct.model(myfitStruct.b_param{i},myfitStruct.Shifts(i),t_mins_t_tips(:,i));
end

%% fit shifted figures;
MDLs = nan(N,length(Shifts),2); %{};
Rsq_mat = nan(N,length(Shifts));
RMSE_cell = cell(length(Shifts),1);
for j = 1:length(Shifts)
    for i = 1:N
        LOGIC_A = t_mins_t_tips(:,i)>t_min & t_mins_t_tips(:,i)<t_max;
        LOGIC = (t_mins_t_tips(:,i)+Shifts(j))>0 & LOGIC_A;
        %--- fit linear
        XX = real(log10(t_mins_t_tips(LOGIC,i)+Shifts(j)));
        YY = real(log10(PhediLocFixed(LOGIC,i)));
        
        [P,~] = polyfit(XX,YY,1);
        %         slp = (YY(end)-YY(1))/(XX(end)-XX(1));
        %         intr = -slp*XX(end)+YY(end);
        %         P = [slp intr];
        
        MDLs(i,j,:) = permute(P(:),[3 2 1]);
        Rsq_mat(i,j) = r_square(YY, polyval(P,XX));
    end
    
end
% mean_Rsq = mean(Rsq_mat(selectedPhedis,:));
mean_Rsq = sum(1-Rsq_mat(selectedPhedis,:));
std_Rsq = std(Rsq_mat(selectedPhedis,:));


%% plot shifted data
FigColors = MyVaryColor(size(PhediLocFixed,2));
n = floor(sqrt(length(Shifts)));
m = ceil(length(Shifts)/n);
% figure;
% ha = tight_subplot(n,m,[.02 .03],[.1 .01],[.01 .01]);

for j = 1:length(Shifts)
    %     axes(ha(j));
    %     hold on;
    figure;
    for i = selectedPhedis
        LOGIC_A = t_mins_t_tips(:,i)>t_min & t_mins_t_tips(:,i)<t_max;
        LOGIC = (t_mins_t_tips(:,i)+Shifts(j))>0 & LOGIC_A;
        hold all;
        %         loglog(t_mins_t_tips(LOGIC,i)+Shifts(j),PhediLocFixed(LOGIC,i),'.-','Color',FigColors(i,:));
        XX = real(log10(t_mins_t_tips(LOGIC,i)+Shifts(j)));
        YY = real(log10(PhediLocFixed(LOGIC,i)));
        %         plot(log10(t_mins_t_tips(LOGIC,i)+Shifts(j)),log10(PhediLocFixed(LOGIC,i)),...
        plot(XX,YY,'.-','Color',FigColors(i,:));
    end
    hold all;
    x = linspace(t_min,t_max,50)+Shifts(j);
    y = sqrt(x);
    Inter_i = mean(MDLs(:,j,2));
    slope_i = mean(MDLs(:,j,1));
    plot(log10(x),log10(.001*y),'k-o','LineWidth',1.5);
    plot(log10(x),Inter_i+slope_i*log10(x),'c-o','LineWidth',1.5);
    title(['shift=',num2str(Shifts(j))]);
    xlabel('shifted time (log10)'); ylabel('location (log10)');
end

%% plot fits statistics
% figure; plot(Shifts,mean_Rsq,'*');
% xlabel('shift'); ylabel('mean RMSE');
% figure; plot(Shifts,std_Rsq,'*');
% xlabel('shift'); ylabel('std RMSE');

% Close2Half = cellfun(@(B) mean(cellfun(@(A) A.Coefficients{2,1},B))-0.5,MDLs);
% figure; plot(Shifts,Close2Half,'*');
% xlabel('shift'); ylabel('power-0.5');

figure; errorbar(Shifts,mean_Rsq,std_Rsq/2,'o');

%% fit parabola and find optimal shift
[~,M] = min(mean_Rsq);
p = polyfit(Shifts(M-2:M+2),mean_Rsq(M-2:M+2),2);
minShift = -p(2)/(2*p(1));
% minShift = find_x0_at_y0(0,Close2Half,Shifts);

%% build model for optimal shift
mdl_final = [];
Rsq_final = [];
for i = 1:N
    %--- plot
    %     LOGIC = t_mins_t_tips(:,j)>t_min & t_mins_t_tips(:,j)<t_max;
    %     %--- fit linear
    %     mdl{j} = fitlm(log10(t_mins_t_tips(LOGIC,j)+minShift),log10(PhediLocFixed(LOGIC,j)));
    %     intercept(j) = mdl{j}.Coefficients{1,1};
    %     x1(j) = mdl{j}.Coefficients{2,1};
    
    LOGIC_A = t_mins_t_tips(:,i)>t_min & t_mins_t_tips(:,i)<t_max;
    LOGIC = (t_mins_t_tips(:,i)+Shifts(j))>0 & LOGIC_A;
    %--- fit linear
    XX = real(log10(t_mins_t_tips(LOGIC,i)+minShift));
    YY = real(log10(PhediLocFixed(LOGIC,i)));
    
    %         [P,~] = polyfit(XX,YY,1);
    slp = (YY(end)-YY(1))/(XX(end)-XX(1));
    intr = -slp*XX(end)-YY(end);
    P = [slp intr];
    
    mdl_final = cat(2,mdl_final,P(:));
    Rsq_final = cat(1,Rsq_final,r_square(YY, polyval(P,XX)));
end
interceptMean = mean(mdl_final(2,selectedPhedis));
x1Mean = mean(mdl_final(1,selectedPhedis));

%% plot fixed data wit sqrt
figure;
xx = linspace(0,max(t_mins_t_tips(:)),500);
yy_sqrt = [zeros(1,499),(10^-interceptMean)*sqrt(xx)];
yy_fit = [zeros(1,499),(10^-interceptMean)*((xx).^x1Mean)];
xx = [-fliplr(xx(2:end)),xx];
% N = size(t_mins_t_tips,2);
% FigColors = MyVaryColor(N);
Markers = repmat({'*'},1,N);
Markers(DataStruct{EvNum}.PhediData.slopeIncline==-1) = {'.'};
Markers(DataStruct{EvNum}.PhediData.slopeIncline==0) = {'diamond'};

for i=selectedPhedis
    hold on;
    plot(t_mins_t_tips(:,i),PhediLocFixed(:,i),...
        '.-','Color',FigColors(i,:),'Marker',Markers{i});
end
hold on;
plot(xx-minShift,yy_fit,'k','LineWidth',1.3);
plot(xx-minShift,yy_sqrt,'Color',rgb('RosyBrown'),'LineWidth',1.3);
