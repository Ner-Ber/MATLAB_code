%% create figure of system
[Cd, Cs, Cr, nu, ro, E, mu, ~,~,~,~]=CrackSolutionMaterialProperties;
%% load for sasmples

% EvNum = 4;  % low Cf
EvNum = 15;  % low Cf
% EvNum = 12;  % high Cf
FontSize = 10;
model = 'exp';
PhediData = PhediStructCell{EvNum}.PhediData;
retrnWave4Cutoff = Cr;

%% slip vel and A/A0
Vx_and_A_F = figure;
AvgPhediStruct = phedi_averagePhedi(PhediStructCell{EvNum});
hold on;
%--plot shaded variance
XX = [AvgPhediStruct.X_mean;flipud(AvgPhediStruct.X_mean)];
YY = [AvgPhediStruct.Vel_mean_x-sqrt(AvgPhediStruct.Vel_var_x);flipud(AvgPhediStruct.Vel_mean_x+sqrt(AvgPhediStruct.Vel_var_x))];
YY = YY(~isnan(XX));    %-- eliminate Nan
XX = XX(~isnan(XX));
XX = XX(~isnan(YY));
YY = YY(~isnan(YY));
% ErrCloud = fill(XX,YY,rgb('PeachPuff'),'linestyle','none','FaceAlpha',0.5,'DisplayName','phedi Vel std');

%--- Add Err Bars
ErrBarLocs = [0.0035 -0.006 -0.02 -0.05];
[~,ErrBarCoor] = min(abs(bsxfun(@minus, AvgPhediStruct.X_mean(:),ErrBarLocs(:)')));
Err = sqrt(AvgPhediStruct.Vel_var_x(ErrBarCoor));
ErrBars = errorbar(AvgPhediStruct.X_mean(ErrBarCoor),AvgPhediStruct.Vel_mean_x(ErrBarCoor),Err,'.','Color',rgb('Khaki'),'LineWidth',2);

%-- select only data before relected waves
[returnLogical, returnIdx] = cutSignalsFromReturningWaves(retrnWave4Cutoff,PhediStructCell{EvNum}, AvgPhediStruct);
[x_Phedi,vel_phedi] = deal(AvgPhediStruct.X_mean(returnLogical),AvgPhediStruct.Vel_mean_x(returnLogical));

%-- plot mean phedi
[~,Ic] = min(abs(Cf_vec(relevantEvents==EvNum)-relevant_colorscale));
vx_color = relevantColorMap(Ic,:);
meanCurve = plot(x_Phedi,vel_phedi,'-','Color',vx_color,'LineWidth',3,'DisplayName','mean phedi');
%-- plot LEFM on interface
LEFMcurve = plot(PhediStructCell{EvNum}.solAtInter.x,PhediStructCell{EvNum}.solAtInter.vx,'k-','LineWidth',2,'DisplayName','LEFM @intrface');

%--- make plot pretty
xlim([-0.05 0.05]);
ylim([-0.02 1.3*max(YY(abs(XX)<=0.05))])
xlabel('x-x_{tip} [m]');
ylabel('velocity [m/s]');
title([PhediStructCell{EvNum}.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf),' vel'])
hold off;


yyaxis right;
%-- superimpose Contact area
x_minXtip_4spatial = PhediStructCell{EvNum}.BigPicRotStruct.x-PhediData.PhotoLocation;
frontTime = PhediStructCell{EvNum}.BigPicRotStruct.frontTime_interp./PhediStructCell{EvNum}.BigPicRotStruct.fps;
[~,I1] = min(abs(x_minXtip_4spatial));
t_minTtip_4tmprl = PhediStructCell{EvNum}.BigPicRotStruct.t-frontTime(I1);
[~,I2] = min(abs(t_minTtip_4tmprl));
A_spatial = PhediStructCell{EvNum}.BigPicRotStruct.DataMatNorm(I2,:);
A_temporal = PhediStructCell{EvNum}.BigPicRotStruct.DataMatNorm(:,I1);
x_minXtip_4tmprl = signal_ChangeTime2Space_mapping(t_minTtip_4tmprl,PhediData.PhotoLocation,PhediStructCell{EvNum}.BigPicRotStruct);

hold on;
plot(x_minXtip_4tmprl,A_temporal,'.-','Color',rgb('DeepSkyBlue'),'DisplayName','pix over time');
% plot(x_minXtip_4spatial,A_spatial,'.-','Color',rgb('DodgerBlue'),'DisplayName','frame over space');
ylabel('A/A_0');
hold off;

xlim([-0.06 0.02]*1e-0);
ylim([0.78 1.05]*1e-0);
Vx_and_A_F.Name = 'Vx_and_A';
Vx_and_A_F.Children.YAxis(2).Color = rgb('DodgerBlue');


%% inset v_peak vs_cf
% A1 = load('G:\Frics\2018-10-10\DAandREsidualsStruct.mat');
% A2 = load('G:\Frics\2018-8-29\DAandREsidualsStruct.mat');
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
%--- extract fracture energy
G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
G_Cell = G_Cell(:);
G_mat = cell2mat(G_Cell);
G_filter = G_mat(Logical_tot);
%--- define display logical
Logical_tot = Cf_logical & damaged_events_logic;


%-- get tau_nl
yieldSxy2ndAll_cell = arrayfun(@(S) S.yieldSxy2ndAll,myStruct,'UniformOutput',0);
yieldSxy2ndAll_cell = yieldSxy2ndAll_cell(:);
yieldSxy2ndAll = cell2mat(yieldSxy2ndAll_cell);



F5b = figure;
F5b.WindowStyle = 'normal';
F5b.Units = 'centimeters';
F5b.Position([3,4]) = [8.6 8.6];
F5b.Name = 'tauVsCf';
tauVsCf = axes;
hold on;
plot(Cf_mat(Logical_tot),yieldSxy2ndAll(Logical_tot),'go','DisplayName','\tau_{nl} (all slopes)');
tauVsCf.YAxis.Scale='log';
axis([[0 1300]*1 10.^[4.9 6.6]]);
xlabel('C_f/C_R');
ylabel('\tau_{nl} [Pa]');

TickVals = [0 0.25 0.5 0.75 1]*Cr;
tauVsCf.XTick= TickVals;
tauVsCf.XTickLabel = cat(2,{0 0.25 0.5 0.75},'C_R');

%--- tau peak by model
Vpeak_All_filter = Vpeak_All(Logical_tot);
tau_p_vec = nan(size(G_filter));
peakVel_vec = nan(size(G_filter));
for i=1:length(G_filter)
    if Cf_filter(i)>=Cr
        thisCf = Cr-1;
    else
        thisCf = Cf_filter(i);
    end
    [X_peak, peakVel_vec(i)] = cohesive_findPeakOfVelByIterations(thisCf/Cr,1e-9,G_filter(i),'L',L*10,model);
    tau_p_vec(i) = CohesiveSolutionGetTau_p(thisCf/Cr,G_filter(i),L*10,model);
end
Ratio = tau_p_vec./peakVel_vec;
tau_p_semiPredictor = Vpeak_All_filter.*Ratio;
slmRatVsCf = cohesive_tau_p_predictorByVx(model);
Ratio_4Cf = ppval(Cf_filter,slmRatVsCf);
tau_p_byModel = Ratio_4Cf.*Vpeak_All_filter;

%-- add data to main plot
plot(Cf_filter,tau_p_byModel,'ro','DisplayName','simple model (all slopes)');

%% inset compare two different taus
F5b_inset = figure; 
F5b_inset.Name = 'comprTau';
compareTau = axes;
hold on
plot(yieldSxy2ndAll(Logical_tot),tau_p_byModel,'o');
plot([0 4e6],[0 4e6],'k-');
ylabel('\tau_p by model');
xlabel('\tau_{nl}');
axis equal;

