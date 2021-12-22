%% create figure of system
[Cd, Cs, Cr, nu, ro, E, mu, ~,~,~,~]=CrackSolutionMaterialProperties;
smthWindow = 65;
maxSrchRegn = 0.03;
Gamma = 3.5;
L = 10*3.2e-3; %3*3.2e-3;
model = 'exp'; %'MySigmoid';
x_vec = -0.03:1e-5:0.01;
retrnWave4Cutoff = Cr;

%% load for sasmples
load('G:\Frics\2018-10-10\allPhediStructures_17186.mat');
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
for iEv = 1:length(relevantEvents)
    PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
    PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
end

FontSize = 14;

displayLowEv = [6, 9, 15, 21];
displayHghEv = [12, 16, 17];
displayCmprEv = [9, 16];
displayCmpr2Mdl = [15, 12];

%% figure for colorbar
CBarFig = figure('Name','CBar');
relevant_colorscale = linspace(0,Cr,1e3);
relevantColorMap = colormap(MyVaryColor(1e3,flipud(My_colorMap)));
colormap(relevantColorMap);
CBar = colorbar;
TickVals = [0 0.25 0.5 0.75 1]*Cr/1260;
CBar.Ticks = TickVals;
CBar.TickLabels = cat(2,{0 0.25 0.5 0.75},'C_R');
CBar.Label.String = 'C_f/C_R ';

%% superimpose mean vx
vx_superImps_F = figure('Name','vx_superImps');
hold on;
Cf_vec = [];
for iEv = 1:length(relevantEvents)
    Cf_vec(iEv) = PhediStructCell{relevantEvents(iEv)}.PhediData.Cf;
end
relevantEvents(Cf_vec>1260) = [];
Cf_vec(Cf_vec>1260) = [];

MeetFrntAt = nan(1,length(relevantEvents));
Cf_save = nan(1,length(relevantEvents));
for iEv = 1:length(relevantEvents)
    Name = [PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details,' Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)];
    AvgPhediStruct = phedi_averagePhedi(PhediStructCell{relevantEvents(iEv)});
    x_Phedi = AvgPhediStruct.X_mean;
    vel_phedi = AvgPhediStruct.Vel_mean_x;
    
    %-- select only data before relected waves
    [returnLogical, returnIdx] = cutSignalsFromReturningWaves(retrnWave4Cutoff,PhediStructCell{relevantEvents(iEv)}, AvgPhediStruct);
    [x_Phedi,vel_phedi] = deal(x_Phedi(returnLogical),vel_phedi(returnLogical));
    
    MeetFrntAt(iEv)= AvgPhediStruct.X_mean(returnIdx);
    Cf_save(iEv) = PhediStructCell{relevantEvents(iEv)}.PhediData.Cf;
    
    
    
    %     WaveTrvlDist = 2*(0.15-PhediStructCell{relevantEvents(iEv)}.PhediData.PhotoLocation);
    %     WaveTrvlTime = WaveTrvlDist/Cr;
    %     [~,returnIdx] = min(abs(AvgPhediStruct.T_mean-WaveTrvlTime));
    %     returnLogical = x_Phedi>=AvgPhediStruct.X_mean(returnIdx);
    %     [x_Phedi_crp,vel_phedi] = deal(x_Phedi(returnLogical),vel_phedi(returnLogical));
    %     x_Phedi = x_Phedi_crp;
    
    
    
    
    
    [~,I] = min(abs(Cf_vec(iEv)-relevant_colorscale));
    plot(x_Phedi,vel_phedi,...
        '-','Color',relevantColorMap(I,:));
    %     plot(x_Phedi,vel_phedi,...
    %             '.-');
end
xlim([-0.08 0.02]*1e-0);
ylim([-0.02 0.55]);
xlabel('x [m]');
ylabel('v_x [m/s]');

%% inset vpeak and cmpr to model
% A1 = load('G:\Frics\2018-10-10\DAandREsidualsStruct.mat');
% A2 = load('G:\Frics\2018-8-29\DAandREsidualsStruct.mat');
A1 = load('E:\Frics\2018-10-10\DAandREsidualsStruct.mat');
A2 = load('E:\Frics\2018-8-29\DAandREsidualsStruct.mat');
myStruct20181010 = A1.myStruct;
myStruct20180829 = A2.myStruct;
myStruct = myStruct20180829;
myStruct((length(myStruct20180829)+1):(length(myStruct20180829)+length(myStruct20181010))) = myStruct20181010(:);
%--- clear empty fields
emptyLogic = arrayfun(@(A) isempty(A.Cf_vec) , myStruct);
myStruct(emptyLogic)=[];

%--- set velocity scale
Cf_scale = linspace(0,1260,9);
CfVals = movmean(Cf_scale,2,'Endpoints','discard');
Colors = MyVaryColor(length(CfVals),flipud(My_colorMap));

%--- create cell of event details
Events_vec = arrayfun(@(S) S.Events,myStruct,'UniformOutput',0);
Events_vec = cell2mat(Events_vec);
Events_vec = Events_vec(:);
Events_cell = num2cell(Events_vec);
Names_cell = arrayfun(@(S) repmat({S.Name},length(S.Events),1),myStruct,'UniformOutput',0);
Names_cell = Names_cell(:);
Names_cell_combined = {};
for i=1:length(Names_cell)
    Names_cell_combined = cat(1,Names_cell_combined,Names_cell{i});
end
event_details = cellfun(@(A,B) [A,' ',num2str(B)],Names_cell_combined,Events_cell,'UniformOutput',0);

%--- define damaged events
damaged_events20180829 = {...           % for 29-8-2018
    '12-8-52 15','12-8-52 17',...
    '14-2-11 6','14-2-11 12',...
    '14-45-14 10','14-45-14 19','14-45-14 28','14-45-14 12',...
    '15-9-25 19',...
    '15-28-23 29','15-28-23 19','15-28-23 18','15-28-23 23','15-28-23 13','15-28-23 8',...
    '18-1-9 14','18-1-9 11',...
    '18-20-0 24','18-20-0 19','18-20-0 7','18-20-0 4'...
    };
damaged_events_wDate20180829 = cellfun(@(A) ['2018-8-29 ',A],damaged_events20180829,'UniformOutput',0);

damaged_events20181010 = {...           % for 10-10-2018
    '12-2-22 28','12-2-22 20',...
    '14-28-36 4',...
    '14-48-1 12','14-48-1 21',...
    '20-28-52 14',...
    '20-48-22 16','20-48-22 3','20-48-22 5',...
    };
damaged_events_wDate20181010 = cellfun(@(A) ['2018-10-10 ',A],damaged_events20181010,'UniformOutput',0);

damaged_events = damaged_events_wDate20180829(:);
damaged_events(length(damaged_events_wDate20180829)+(1:length(damaged_events_wDate20181010))) = damaged_events_wDate20181010(:);
damaged_events_logic = ~ismember(event_details,damaged_events);

N = cellfun(@(S) min(strfind(S,' ')), event_details,'UniformOutput',0);
event_dates = cellfun(@(A,n) A(1:(n-1)),event_details,N,'UniformOutput',0);
[~,~,groupCatagrs0] = unique(event_dates);
presentedData = cellfun(@(A) ~isempty(strfind(A,'17-18-6')),Names_cell_combined);
groupCatagrs = groupCatagrs0;
groupCatagrs(presentedData) = 3;
%--- extract Cf
Cf_Cell = arrayfun(@(S) S.Cf_vec,myStruct,'UniformOutput',0);
Cf_Cell = Cf_Cell(:);
Cf_mat = cell2mat(Cf_Cell);
Cf_logical = Cf_mat<=1260 & Cf_mat>=0;

%--- define display logical
Logical_tot = Cf_logical & damaged_events_logic;

%--- exctract peak velocity
Vpeak_All_cell = arrayfun(@(S) S.peakVel_vecAll,myStruct,'UniformOutput',0);
Vpeak_All_cell = Vpeak_All_cell(:);
Vpeak_All = cell2mat(Vpeak_All_cell);

%%-- peak velocity by model
Gamma = 3.5;
Cf_vec_4Model = linspace(1^2.8,1250^2.8,100).^(1/2.8);
L = 10*3.2e-3;
x_vec = [-2e-3:1e-4:(-8e-4-5e-6),-8e-4:5e-6:2e-4]; %-0.03:1e-5:0.01;
L_estimate_model = CohesiveSolutionGetL(Cf_vec_4Model./Cr,Gamma,2e6)/10;
maxCohesive = nan(size(Cf_vec_4Model));
for c = 1:length(Cf_vec_4Model)
    %     sol_cohesive = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_vec_4Model(c)/Cr,1e-9,Gamma,'L',L,model,[],x_vec);
    %     [maxCohesive(c),~] = max(sol_cohesive.vx(:));
    [~, maxCohesive(c)] = cohesive_findPeakOfVelByIterations(Cf_vec_4Model(c)/Cr,1e-9,Gamma,'L',L,model);
end

%--- prepare and plot
Vpeak_All_filter = Vpeak_All(Logical_tot);
Cf_filter = Cf_mat(Logical_tot);
groupCatagrs_filter = groupCatagrs(Logical_tot);

% axes(v_peak_vs_Cf); hold on;
v_peak_vs_Cf_F = figure('Name','v_peak_vs_Cf');
% v_peak_vs_Cf = v_peak_vs_Cf_F.Children;
hold on;
plot(Cf_vec_4Model,maxCohesive,'Color',[1 0.8 0.2], 'LineWidth', 1.5);
plot(Cf_filter(groupCatagrs_filter==2),Vpeak_All_filter(groupCatagrs_filter==2),'ks','LineWidth',1,'MarkerSize',4);
plot(Cf_filter(groupCatagrs_filter==1),Vpeak_All_filter(groupCatagrs_filter==1),'ko','LineWidth',1,'MarkerSize',4);
plot(Cf_filter(groupCatagrs_filter==3),Vpeak_All_filter(groupCatagrs_filter==3),'ko','MarkerFaceColor','r','LineWidth',1,'MarkerSize',4);
xlim([0 1300]);
ylim([0 1.4]);
xlabel('C_f/C_R');
ylabel('v_x peak [m/s]');
v_peak_vs_Cf_F.Children.XAxis.TickValues = [0:0.5:1]*Cr;
v_peak_vs_Cf_F.Children.XAxis.TickLabels= cat(2,num2cell([0 0.5]),'C_R');

%% set Cf scale and colors
% Cf_vec = [];
% for iEv = 1:length(relevantEvents)
%     Cf_vec(iEv) = PhediStructCell{relevantEvents(iEv)}.PhediData.Cf;
% end
% relevantEvents(Cf_vec>1260) = [];
% Cf_vec(Cf_vec>1260) = [];
% relevant_colorscale = linspace(0,Cr,1e3);
% relevantColorMap = colormap(MyVaryColor(1e3,flipud(My_colorMap)));

%% plot for low regime
vx_lowCfClps_F = figure('Name','lowCfCollapse');
% vx_lowCfClps = axes;
hold on;
for iEv = find(ismember(relevantEvents,relevantEvents))
    AvgPhediStruct = phedi_averagePhedi(PhediStructCell{relevantEvents(iEv)});
    x_Phedi_all = AvgPhediStruct.X_mean;
    vel_phedi_all = AvgPhediStruct.Vel_mean_x;
    vel_phedi = vel_phedi_all(abs(x_Phedi_all)<=0.06);
    x_Phedi = x_Phedi_all(abs(x_Phedi_all)<=0.06);
    
    [~,I] = min(abs(Cf_vec(iEv)-relevant_colorscale));
    
    %-- smooth and normalize
    Y = smooth(vel_phedi,smthWindow);
    vel_phedi_norm = Y/max(Y(abs(x_Phedi)<=maxSrchRegn));
    
    %-- select only data before relected waves
    
    [returnLogical, returnIdx] = cutSignalsFromReturningWaves(retrnWave4Cutoff,PhediStructCell{relevantEvents(iEv)}, AvgPhediStruct);
    returnLogical = returnLogical(abs(x_Phedi_all)<=0.06);
    [x_Phedi,vel_phedi_norm] = deal(x_Phedi(returnLogical),vel_phedi_norm(returnLogical));
    
    plot(x_Phedi,vel_phedi_norm,...
        '-','Color',relevantColorMap(I,:),...
        'DisplayName',[PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details, '  Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)]);
end
colormap(relevantColorMap);
xlim([-0.055 0.055]*1e-0);
ylim([-0.05 1.3]*1e-0);
xlabel('x-x_{tip} [m]');
ylabel('v_{x, normalized} [m/s]');

%% plot for high regime
figure('Name','highCfCollapse');
vx_highCfClps = axes;
hold on;
for iEv = find(ismember(relevantEvents,displayHghEv))
    AvgPhediStruct = phedi_averagePhedi(PhediStructCell{relevantEvents(iEv)});
    x_Phedi_all = AvgPhediStruct.X_mean;
    vel_phedi_all = AvgPhediStruct.Vel_mean_x;
    vel_phedi = vel_phedi_all(abs(x_Phedi_all)<=0.06);
    x_Phedi = x_Phedi_all(abs(x_Phedi_all)<=0.06);
    
    
    [~,I] = min(abs(Cf_vec(iEv)-relevant_colorscale));
    
    %-- smooth and normalize
    Y = smooth(vel_phedi,smthWindow);
    vel_phedi_norm = Y/max(Y(abs(x_Phedi)<=maxSrchRegn));
    
    %-- select only data before relected waves
    [returnLogical, returnIdx] = cutSignalsFromReturningWaves(retrnWave4Cutoff,PhediStructCell{relevantEvents(iEv)}, AvgPhediStruct);
    returnLogical = returnLogical(abs(x_Phedi_all)<=0.06);
    [x_Phedi,vel_phedi_norm] = deal(x_Phedi(returnLogical),vel_phedi_norm(returnLogical));
    
    plot(x_Phedi,vel_phedi_norm,...
        '-','Color',relevantColorMap(I,:),...
        'DisplayName',[PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details, '  Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)]);
end
colormap(relevantColorMap);
xlim([-0.055 0.055]*1e-0);
ylim([-0.05 1.3]*1e-0);
xlabel('x-x_{tip} [m]');
%% plot for comparison sublot
figure('Name','CmprVelCollapse');
vx_CmprCf = axes;
hold on;
for iEv = find(ismember(relevantEvents,displayCmprEv))
    AvgPhediStruct = phedi_averagePhedi(PhediStructCell{relevantEvents(iEv)});
    x_Phedi_all = AvgPhediStruct.X_mean;
    vel_phedi_all = AvgPhediStruct.Vel_mean_x;
    vel_phedi = vel_phedi_all(abs(x_Phedi_all)<=0.06);
    x_Phedi = x_Phedi_all(abs(x_Phedi_all)<=0.06);
    
    [~,I] = min(abs(Cf_vec(iEv)-relevant_colorscale));
    
    %-- smooth and normalize
    Y = smooth(vel_phedi,smthWindow);
    vel_phedi_norm = Y/max(Y(abs(x_Phedi)<=maxSrchRegn));
    
    %-- select only data before relected waves
    [returnLogical, returnIdx] = cutSignalsFromReturningWaves(retrnWave4Cutoff,PhediStructCell{relevantEvents(iEv)}, AvgPhediStruct);
    returnLogical = returnLogical(abs(x_Phedi_all)<=0.06);
    [x_Phedi,vel_phedi_norm] = deal(x_Phedi(returnLogical),vel_phedi_norm(returnLogical));
    
    plot(x_Phedi,vel_phedi_norm,...
        '-','Color',relevantColorMap(I,:),...
        'DisplayName',[PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details, '  Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)]);
end
% CBar = colorbar;
% colormap(relevantColorMap);
xlim([-0.055 0.055]*1e-0);
ylim([-0.05 1.3]*1e-0);
xlabel('x-x_{tip} [m]');
%% plot for inset - compare to model
x_vec = -0.06:1e-4:0.02;
figure('Name','Cmpr2Model');
Cmpr2Mdl = axes; hold on;
for iEv = find(ismember(relevantEvents,displayCmpr2Mdl))
    % for iEv = 1:length(Cf_vec)
    AvgPhediStruct = phedi_averagePhedi(PhediStructCell{relevantEvents(iEv)});
    x_Phedi = AvgPhediStruct.X_mean;
    vel_phedi = AvgPhediStruct.Vel_mean_x;
    [~,I] = min(abs(Cf_vec(iEv)-relevant_colorscale));
    
    %-- smooth and normalize
    Y = smooth(vel_phedi,smthWindow);
    vel_phedi_norm = Y/max(Y(abs(x_Phedi)<=maxSrchRegn));
    
    %-- select only data before relected waves
    [returnLogical, returnIdx] = cutSignalsFromReturningWaves(retrnWave4Cutoff,PhediStructCell{relevantEvents(iEv)}, AvgPhediStruct);
%     [x_Phedi,vel_phedi] = deal(x_Phedi(returnLogical),vel_phedi(returnLogical));
    [x_Phedi,vel_phedi] = deal(x_Phedi(returnLogical),Y(returnLogical));
    
    plot(x_Phedi,vel_phedi,...
        '.-','Color',relevantColorMap(I,:),...
        'DisplayName',[PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details, '  Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)]);
    
    %     sol_cohesive = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_vec(iEv)/Cr,1e-9,ppval(GammaFit,Cf_vec(iEv)),'L',L,model,[],x_vec);
    sol_cohesive = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_vec(iEv)/Cr,1e-9,Gamma,'L',L,model,[],x_vec);
    plot(sol_cohesive.x,sol_cohesive.vx,'--','LineWidth',2,'Color',relevantColorMap(I,:));
end
Cmpr2Mdl.Children(1).Marker = '.';
Cmpr2Mdl.Children(3).Marker = '.';
Cmpr2Mdl.Children(2).LineWidth = 1;
Cmpr2Mdl.Children(4).LineWidth = 1;
Cmpr2Mdl.FontSize = 10;
Cmpr2Mdl.YAxis.TickValues = [0 0.2 0.4];
Cmpr2Mdl.XAxis.TickValues = [-0.05 -0.02 0 0.02];
xlim([-0.08 0.02]*1e-0);
ylim([-0.02 0.55]);
xlabel('x-x_{tip} [m]');
ylabel('v_x [m/s]');

