%% set data
DataStruct = {PhediStruct_3_11};
EvNum = 1;
preKink = -4.6e-4;
myfitStruct = phedi_fitSineToPhediTraj(DataStruct{EvNum}.PhediData,[],2,1);
%-- parameters:
t_min = 0.003/DataStruct{EvNum}.BigPicRotStruct.AvgCf_scratches;
t_max = 0.025/DataStruct{EvNum}.BigPicRotStruct.AvgCf_scratches;
Shifts = linspace(-2*1e-5,2*1e-5,200);
selectedPhedis = [4 6 8 12 14 16 23]; %[2:2:16,19,21,23,25];

%% create fixed matrix
PhediLoc = DataStruct{EvNum}.PhediData.PhediLocation;
t_mins_t_tips = bsxfun(@minus,DataStruct{EvNum}.PhediData.timeVec,DataStruct{EvNum}.PhediData.t_tips(:)');

N = size(t_mins_t_tips,2);
PhediLocFixed = nan(size(PhediLoc));
for i=1:N
    PhediLocFixed(:,i) = PhediLoc(:,i)-myfitStruct.model(myfitStruct.b_param{i},myfitStruct.Shifts(i),t_mins_t_tips(:,i));
end

%% fit shifted figures;

%--- interpolate for common time steps
t_unique = unique(t_mins_t_tips);
T_phedis = repmat(t_unique(:),1,N);
PhediNumbering = repmat(1:N,length(t_unique),1);
% PhediLocFixed_iterped = interp2(t_mins_t_tips,repmat(1:N,size(t_mins_t_tips,1),1),PhediLocFixed,T_phedis,PhediNumbering);
PhediLocFixed_iterped = nan(size(T_phedis));
for i=1:N
    PhediLocFixed_iterped(:,i) = interp1(t_mins_t_tips(:,i),PhediLocFixed(:,i),T_phedis(:,i));
end

%--- create mean measurement
t_mins_t_tips_mean = mean(T_phedis(:,selectedPhedis),2);
PhediLocFixed_mean = mean(PhediLocFixed_iterped(:,selectedPhedis),2);

%--- create model for mean measurement
MDLs = nan(1,length(Shifts),2); %{};
Rsq_mat = nan(1,length(Shifts));
for j = 1:length(Shifts)
    LOGIC_A = t_mins_t_tips_mean(:)>t_min & t_mins_t_tips_mean(:)<t_max;
    LOGIC = (t_mins_t_tips_mean(:)+Shifts(j))>0 & LOGIC_A;
    %--- fit linear
    XX = real(log10(t_mins_t_tips_mean(LOGIC)+Shifts(j)));
    YY = real(log10(PhediLocFixed_mean(LOGIC)));
    
    [P,~] = polyfit(XX,YY,1);
    %         slp = (YY(end)-YY(1))/(XX(end)-XX(1));
    %         intr = -slp*XX(end)+YY(end);
    %         P = [slp intr];
    
    MDLs(1,j,:) = permute(P(:),[3 2 1]);
    Rsq_mat(j) = r_square(YY, polyval(P,XX));
end
% mean_Rsq = mean(Rsq_mat(selectedPhedis,:));
mean_Rsq = 1-Rsq_mat;


%% plot shifted data
% FigColors = MyVaryColor(size(PhediLocFixed,2));
% n = floor(sqrt(length(Shifts)));
% m = ceil(length(Shifts)/n);
% % figure;
% % ha = tight_subplot(n,m,[.02 .03],[.1 .01],[.01 .01]);
% 
% for j = 1:length(Shifts)
%     %     axes(ha(j));
%     %     hold on;
%     figure;
%         LOGIC_A = t_mins_t_tips_mean(:)>t_min & t_mins_t_tips_mean(:)<t_max;
%         LOGIC = (t_mins_t_tips_mean(:)+Shifts(j))>0 & LOGIC_A;
%         hold all;
%         %         loglog(t_mins_t_tips(LOGIC,i)+Shifts(j),PhediLocFixed(LOGIC,i),'.-','Color',FigColors(i,:));
%         XX = real(log10(t_mins_t_tips_mean(LOGIC)+Shifts(j)));
%         YY = real(log10(PhediLocFixed_mean(LOGIC)));
%         %         plot(log10(t_mins_t_tips(LOGIC,i)+Shifts(j)),log10(PhediLocFixed(LOGIC,i)),...
%         plot(XX,YY,'.-','Color',FigColors(i,:));
%     hold all;
%     x = linspace(t_min,t_max,50)+Shifts(j);
%     y = sqrt(x);
%     Inter_i = mean(MDLs(:,j,2));
%     slope_i = mean(MDLs(:,j,1));
%     plot(log10(x),log10(.001*y),'k-o','LineWidth',1.5);
%     plot(log10(x),Inter_i+slope_i*log10(x),'c-o','LineWidth',1.5);
%     title(['shift=',num2str(Shifts(j))]);
%     xlabel('shifted time (log10)'); ylabel('location (log10)');
% end
% 
% %% plot fits statistics
% % figure; plot(Shifts,mean_Rsq,'*');
% % xlabel('shift'); ylabel('mean RMSE');
% % figure; plot(Shifts,std_Rsq,'*');
% % xlabel('shift'); ylabel('std RMSE');
% 
% % Close2Half = cellfun(@(B) mean(cellfun(@(A) A.Coefficients{2,1},B))-0.5,MDLs);
% % figure; plot(Shifts,Close2Half,'*');
% % xlabel('shift'); ylabel('power-0.5');
% 
% figure; plot(Shifts,mean_Rsq,'o');
% 
% %% fit parabola and find optimal shift
% [~,M] = min(mean_Rsq);
% p = polyfit(Shifts(M-2:M+2),mean_Rsq(M-2:M+2),2);
% minShift = -p(2)/(2*p(1));
% % minShift = find_x0_at_y0(0,Close2Half,Shifts);
% 
% %% build model for optimal shift
% LOGIC_A = t_mins_t_tips_mean(:)>t_min & t_mins_t_tips_mean(:)<t_max;
% LOGIC = (t_mins_t_tips_mean(:)+minShift)>0 & LOGIC_A;
% %--- fit linear
% XX = real(log10(t_mins_t_tips_mean(LOGIC)+minShift));
% YY = real(log10(PhediLocFixed_mean(LOGIC)));
% 
% [MDL_final,~] = polyfit(XX,YY,1);
% %         slp = (YY(end)-YY(1))/(XX(end)-XX(1));
% %         intr = -slp*XX(end)+YY(end);
% %         P = [slp intr];
% 
% Rsq_mat(j) = r_square(YY, polyval(MDL_final,XX));
% interceptMean = MDL_final(2);
% x1Mean = MDL_final(1);
% 
% %% plot fixed data wit sqrt
% figure;
% xx = linspace(0,max(t_mins_t_tips_mean(:)),500);
% yy_sqrt = [zeros(1,499),(10^interceptMean)*sqrt(xx)];
% yy_fit = [zeros(1,499),(10^interceptMean)*((xx).^x1Mean)];
% xx = [-fliplr(xx(2:end)),xx];
% % N = size(t_mins_t_tips,2);
% % FigColors = MyVaryColor(N);
% Markers = repmat({'*'},1,N);
% Markers(DataStruct{EvNum}.PhediData.slopeIncline==-1) = {'.'};
% Markers(DataStruct{EvNum}.PhediData.slopeIncline==0) = {'diamond'};
% 
% for i=selectedPhedis
%     hold on;
%     plot(t_mins_t_tips(:,i),PhediLocFixed(:,i),...
%         '.-','Color',FigColors(i,:),'Marker',Markers{i});
% end
% hold on;
% plot(xx-minShift,yy_fit,'k','LineWidth',1.3);
% plot(xx-minShift,yy_sqrt,'Color',rgb('RosyBrown'),'LineWidth',1.3);
