%% choose locations
Aloc = -0.04;
VresLoc = -0.05;
[Cd, Cs, Cr, nu, ro, E, mu, ~, PlaneStrain, tau_p, Xc0] = CrackSolutionMaterialProperties;
%% load structure

% A1 = load('F:\NeriFricsData\Frics\2018-10-10\DAandREsidualsStruct.mat');
% A2 = load('F:\NeriFricsData\Frics\2018-8-29\DAandREsidualsStruct.mat');
A1 = load('E:\Frics\2018-10-10\DAandREsidualsStruct.mat');
A2 = load('E:\Frics\2018-8-29\DAandREsidualsStruct.mat');
myStruct20181010 = A1.myStruct;
myStruct20180829 = A2.myStruct;
myStruct = myStruct20180829;
myStruct((length(myStruct20180829)+1):(length(myStruct20180829)+length(myStruct20181010))) = myStruct20181010(:);
%% clear empty fields
emptyLogic = arrayfun(@(A) isempty(A.Cf_vec) , myStruct);
myStruct(emptyLogic)=[];

%% set velocity scale
Cf_scale = linspace(0,1260,9);
CfVals = movmean(Cf_scale,2,'Endpoints','discard');
Colors = MyVaryColor(length(CfVals),flipud(My_colorMap));

%% create cell of event details
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

%% define damaged events
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
[~,~,groupCatagrs] = unique(event_dates);
%% extract Cf
Cf_Cell = arrayfun(@(S) S.Cf_vec,myStruct,'UniformOutput',0);
Cf_Cell = Cf_Cell(:);
Cf_mat = cell2mat(Cf_Cell);
Cf_logical = Cf_mat<=1260 & Cf_mat>=0;

%% define display logical
Logical_tot = Cf_logical & damaged_events_logic;
Cf_filter = Cf_mat(Logical_tot);

%% get data
yieldSxy2ndAll_cell = arrayfun(@(S) S.yieldSxy2ndAll,myStruct,'UniformOutput',0);
yieldSxy2ndAll_cell = yieldSxy2ndAll_cell(:);
yieldSxy2ndAll = cell2mat(yieldSxy2ndAll_cell);

%--- v peak by model
Vpeak_All_cell = arrayfun(@(S) S.peakVel_vecAll,myStruct,'UniformOutput',0);
Vpeak_All_cell = Vpeak_All_cell(:);
Vpeak_All = cell2mat(Vpeak_All_cell);

Vpeak_vec = Vpeak_All(Logical_tot);
tau_nl_vec = yieldSxy2ndAll(Logical_tot);

[Vpk_vec_smth, Tnl_vec_var, Tnl_vec_mean] = my_movVar(Vpeak_vec, tau_nl_vec, 0.3, 0.04);

%--- devide into above and below Cc
Cc = 0.8*Cr;
final_Cf = Cf_mat(Logical_tot);
Cc_logic_below = final_Cf<Cc;
Cc_logic_above = final_Cf>=Cc;
Vpeak_vec_below = Vpeak_vec(Cc_logic_below);
Vpeak_vec_above = Vpeak_vec(Cc_logic_above);
tau_nl_vec_below = tau_nl_vec(Cc_logic_below);
tau_nl_vec_above = tau_nl_vec(Cc_logic_above);

%% add tau_p from model
%--- extract fracture energy
G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
G_Cell = G_Cell(:);
G_mat = cell2mat(G_Cell);
G_filter = G_mat(Logical_tot);
%--- tau peak by model
model = 'exp';
L = 3e-3*10;
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

tau_p_byModel_below = tau_p_byModel(Cc_logic_below);
tau_p_byModel_above = tau_p_byModel(Cc_logic_above);

%% plot
figure;
tauVsCf = axes;
hold on;
% plot(Vpeak_vec,tau_nl_vec*1e-6,'o');
plot(Vpeak_vec_below,tau_nl_vec_below*1e-6,'bo');
plot(Vpeak_vec_above,tau_nl_vec_above*1e-6,'ro');
plot(Vpk_vec_smth,Tnl_vec_mean*1e-6,'k');
%-- add data to main plot
% plot(Vpeak_vec,tau_p_byModel,'ro','DisplayName','simple model (all slopes)');
plot(Vpeak_vec_below,tau_p_byModel_below*1e-6,'.','DisplayName','simple model (all slopes)', 'Color',rgb('SlateGray'), 'MarkerSize',15);
plot(Vpeak_vec_above,tau_p_byModel_above*1e-6,'.','DisplayName','simple model (all slopes)', 'Color',rgb('Peru'), 'MarkerSize',15);

% tauVsCf.XAxis.Scale='log';
% axis([[0 1300]*1 10.^[4.9 6.6]]);
xlabel('v_{peak} [m/s]');
ylabel('\tau_{nl} [MPa]');
pbaspect([5,3,1])
xlim([1e-2,10^(0.5)])
pbaspect([8,5,1])
