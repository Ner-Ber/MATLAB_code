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

% %% txtract Ar
% [~,locLogicA] = min(abs(myStruct(1).residual_locs - Aloc));
% 
% Ar_cell_spt= arrayfun(@(S) S.Ar_vec_spt(locLogicA,:),myStruct,'UniformOutput',0);
% Ar_vec_spt = cell2mat(Ar_cell_spt);
% Ar_cell_tmp= arrayfun(@(S) S.Ar_vec_tmp(locLogicA,:),myStruct,'UniformOutput',0);
% Ar_vec_tmp = cell2mat(Ar_cell_tmp);
% Ar_cell_sml = arrayfun(@(S) S.Ar_vec_sml(locLogicA,:),myStruct,'UniformOutput',0);
% Ar_vec_sml = cell2mat(Ar_cell_sml);
% 
% DA_cell_spt= arrayfun(@(S) S.DA_vec_spt(locLogicA,:),myStruct,'UniformOutput',0);
% DA_vec_spt = cell2mat(DA_cell_spt);
% 
% %% exctract peak velocity
% Vpeak_All_cell = arrayfun(@(S) S.peakVel_vecAll,myStruct,'UniformOutput',0);
% Vpeak_All_cell = Vpeak_All_cell(:);
% Vpeak_All = cell2mat(Vpeak_All_cell);
% Vpeak_Neg_cell = arrayfun(@(S) S.peakVel_vecNeg,myStruct,'UniformOutput',0);
% Vpeak_Neg_cell = Vpeak_Neg_cell(:);
% Vpeak_Neg = cell2mat(Vpeak_Neg_cell);
% Vpeak_Pos_cell = arrayfun(@(S) S.peakVel_vecPos,myStruct,'UniformOutput',0);
% Vpeak_Pos_cell = Vpeak_Pos_cell(:);
% Vpeak_Pos = cell2mat(Vpeak_Pos_cell);
% 
% %% extract residual velocity
% [~,locLogicV] = min(abs(myStruct(1).residual_locs - VresLoc));
% 
% 
% Vresiduals_All_cell = arrayfun(@(S) S.Vresiduals_matAll(locLogicV,:),myStruct,'UniformOutput',0);
% V_res_all = cell2mat(Vresiduals_All_cell);
% Vresiduals_Neg_cell = arrayfun(@(S) S.Vresiduals_matNeg(locLogicV,:),myStruct,'UniformOutput',0);
% V_res_neg = cell2mat(Vresiduals_Neg_cell);
% Vresiduals_Pos_cell = arrayfun(@(S) S.Vresiduals_matPos(locLogicV,:),myStruct,'UniformOutput',0);
% V_res_pos = cell2mat(Vresiduals_Pos_cell);

%% plot
rise2XAllcell = arrayfun(@(S) S.rise2XAll,myStruct,'UniformOutput',0);
rise2XAllcell = rise2XAllcell(:);
rise2XAll_vec = cell2mat(rise2XAllcell);
rise2XNegcell = arrayfun(@(S) S.rise2XNeg,myStruct,'UniformOutput',0);
rise2XNegcell = rise2XNegcell(:);
rise2XNeg_vec = cell2mat(rise2XNegcell);
rise2XPoscell = arrayfun(@(S) S.rise2XPos,myStruct,'UniformOutput',0);
rise2XPoscell = rise2XPoscell(:);
rise2XPos_vec = cell2mat(rise2XPoscell);

figure; ax = axes; hold on;
colormap(flipud(My_colorMap));
C_vec = Cf_mat(Logical_tot)/Cr;
Xnl_vec = rise2XAll_vec(Logical_tot)*1e3;
scatter(C_vec,Xnl_vec,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
% scatter(Cf_mat(Logical_tot)/Cr,rise2XNeg_vec(Logical_tot)*1e3,[],groupCatagrs(Logical_tot),'filled','d','DisplayName','neg slopes');
% scatter(Cf_mat(Logical_tot)/Cr,rise2XPos_vec(Logical_tot)*1e3,[],groupCatagrs(Logical_tot),'filled','s','DisplayName','pos slopes');
[C_vec_smth, Xnl_vec_var, Xnl_vec_mean] = my_movVar(C_vec, Xnl_vec, 0.18, 0.06);
plot(C_vec_smth,Xnl_vec_mean, 'DisplayName','x_{nl} mean');
title('x_{nl} vs. C_f');
ylabel('x_{nl} [mm]');
xlabel('C_f/C_R [m/s]');
pbaspect([5,3,1])