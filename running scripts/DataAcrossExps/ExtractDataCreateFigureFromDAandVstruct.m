%% create extract data and figure from structure
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

%% txtract Ar
[~,locLogicA] = min(abs(myStruct(1).residual_locs - Aloc));

Ar_cell_spt= arrayfun(@(S) S.Ar_vec_spt(locLogicA,:),myStruct,'UniformOutput',0);
Ar_vec_spt = cell2mat(Ar_cell_spt);
Ar_cell_tmp= arrayfun(@(S) S.Ar_vec_tmp(locLogicA,:),myStruct,'UniformOutput',0);
Ar_vec_tmp = cell2mat(Ar_cell_tmp);
Ar_cell_sml = arrayfun(@(S) S.Ar_vec_sml(locLogicA,:),myStruct,'UniformOutput',0);
Ar_vec_sml = cell2mat(Ar_cell_sml);

DA_cell_spt= arrayfun(@(S) S.DA_vec_spt(locLogicA,:),myStruct,'UniformOutput',0);
DA_vec_spt = cell2mat(DA_cell_spt);

%% exctract peak velocity
Vpeak_All_cell = arrayfun(@(S) S.peakVel_vecAll,myStruct,'UniformOutput',0);
Vpeak_All_cell = Vpeak_All_cell(:);
Vpeak_All = cell2mat(Vpeak_All_cell);
Vpeak_Neg_cell = arrayfun(@(S) S.peakVel_vecNeg,myStruct,'UniformOutput',0);
Vpeak_Neg_cell = Vpeak_Neg_cell(:);
Vpeak_Neg = cell2mat(Vpeak_Neg_cell);
Vpeak_Pos_cell = arrayfun(@(S) S.peakVel_vecPos,myStruct,'UniformOutput',0);
Vpeak_Pos_cell = Vpeak_Pos_cell(:);
Vpeak_Pos = cell2mat(Vpeak_Pos_cell);

%% extract residual velocity
[~,locLogicV] = min(abs(myStruct(1).residual_locs - VresLoc));


Vresiduals_All_cell = arrayfun(@(S) S.Vresiduals_matAll(locLogicV,:),myStruct,'UniformOutput',0);
V_res_all = cell2mat(Vresiduals_All_cell);
Vresiduals_Neg_cell = arrayfun(@(S) S.Vresiduals_matNeg(locLogicV,:),myStruct,'UniformOutput',0);
V_res_neg = cell2mat(Vresiduals_Neg_cell);
Vresiduals_Pos_cell = arrayfun(@(S) S.Vresiduals_matPos(locLogicV,:),myStruct,'UniformOutput',0);
V_res_pos = cell2mat(Vresiduals_Pos_cell);

%% plot Ar vs. Vres
if 0
    colormap(flipud(My_colorMap));
    figure; ax = axes; hold on;
    % semilogx(V_res_all(Logical_tot),Ar_vec_spt(Logical_tot),'*','DisplayName','spatial, all slopes');
    % semilogx(V_res_all(Logical_tot),Ar_vec_tmp(Logical_tot),'*','DisplayName','temporal, all slopes');
    % semilogx(V_res_all(Logical_tot),Ar_vec_sml(Logical_tot),'*','DisplayName','small, all slopes');
    % semilogx(V_res_neg(Logical_tot),Ar_vec_spt(Logical_tot),'s','DisplayName','spatial, neg slopes')
    % semilogx(V_res_neg(Logical_tot),Ar_vec_tmp(Logical_tot),'s','DisplayName','temporal, neg slopes');
    % semilogx(V_res_neg(Logical_tot),Ar_vec_sml(Logical_tot),'s','DisplayName','small, neg slopes');
    % semilogx(V_res_pos(Logical_tot),Ar_vec_spt(Logical_tot),'v','DisplayName','spatial, pos slopes')
    % semilogx(V_res_pos(Logical_tot),Ar_vec_tmp(Logical_tot),'v','DisplayName','temporal, pos slopes');
    % semilogx(V_res_pos(Logical_tot),Ar_vec_sml(Logical_tot),'v','DisplayName','small, pos slopes');
    
    scatter(V_res_all(Logical_tot),Ar_vec_spt(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','spatial, all slopes');
    scatter(V_res_all(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','temporal, all slopes');
    scatter(V_res_all(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','small, all slopes');
    scatter(V_res_neg(Logical_tot),Ar_vec_spt(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','spatial, neg slopes')
    scatter(V_res_neg(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','temporal, neg slopes');
    scatter(V_res_neg(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','small, neg slopes');
    scatter(V_res_pos(Logical_tot),Ar_vec_spt(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','spatial, pos slopes')
    scatter(V_res_pos(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','temporal, pos slopes');
    scatter(V_res_pos(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','o','MarkerEdgeColor','k','DisplayName','small, pos slopes');
    
    ylabel('A_r/A_0');
    xlabel('v_{res} [m/s]');
    ax.XAxis.Scale='log';
    title('A_r vs. v_{res}, res location=-0.03');
end

%% plot Ar vs. Vpeak
if 0
    colormap(flipud(My_colorMap));
    figure; ax = axes; hold on;
    
    %     scatter(Vpeak_All(Logical_tot),Ar_vec_spt(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','spatial, all slopes');
    %     scatter(Vpeak_All(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','temporal, all slopes');
    %     scatter(Vpeak_All(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','small, all slopes');
    %     scatter(Vpeak_Neg(Logical_tot),Ar_vec_spt(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','spatial, neg slopes')
    %     scatter(Vpeak_Neg(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','temporal, neg slopes');
    %     scatter(Vpeak_Neg(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','small, neg slopes');
    %     scatter(Vpeak_Pos(Logical_tot),Ar_vec_spt(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','spatial, pos slopes')
    %     scatter(Vpeak_Pos(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','temporal, pos slopes');
    %     scatter(Vpeak_Pos(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','small, pos slopes');
    
    scatter(Cf_mat(Logical_tot),Ar_vec_spt(Logical_tot),[],Vpeak_All(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','spatial, all slopes');
    %     scatter(Cf_mat(Logical_tot),DA_vec_spt(Logical_tot),[],Vpeak_All(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','spatial, all slopes');
    %     scatter(Vpeak_All(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','temporal, all slopes');
    %     scatter(Vpeak_All(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','small, all slopes');
    %     scatter(Vpeak_Neg(Logical_tot),Ar_vec_spt(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','spatial, neg slopes')
    %     scatter(Vpeak_Neg(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','temporal, neg slopes');
    %     scatter(Vpeak_Neg(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','small, neg slopes');
    %     scatter(Vpeak_Pos(Logical_tot),Ar_vec_spt(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','spatial, pos slopes')
    %     scatter(Vpeak_Pos(Logical_tot),Ar_vec_tmp(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','temporal, pos slopes');
    %     scatter(Vpeak_Pos(Logical_tot),Ar_vec_sml(Logical_tot),[],Cf_mat(Logical_tot),'filled','Marker','s','MarkerEdgeColor','k','DisplayName','small, pos slopes');
    
    %     ylabel('A_r/A_0');
    ylabel('\DeltaA/A_0');
    %     xlabel('v_{peak} [m/s]');
    xlabel('C_{f} [m/s]');
    ax.XAxis.Scale='log';
    %     title(['A_r vs. v_{peak}, res loc=',num2str(VresLoc),' A loc=',num2str(Aloc)]);
    title(['\DeltaA vs. C_{f}, res loc=',num2str(VresLoc),' A loc=',num2str(Aloc)]);
end

%% plot FWHM
if 0
    FWHMall_cell = arrayfun(@(S) S.FWHMall,myStruct,'UniformOutput',0);
    FWHMall_cell = FWHMall_cell(:);
    FWHMall_vec = cell2mat(FWHMall_cell);
    FWHMneg_cell = arrayfun(@(S) S.FWHMneg,myStruct,'UniformOutput',0);
    FWHMneg_cell = FWHMneg_cell(:);
    FWHMneg_vec = cell2mat(FWHMneg_cell);
    FWHMpos_cell = arrayfun(@(S) S.FWHMpos,myStruct,'UniformOutput',0);
    FWHMpos_cell = FWHMpos_cell(:);
    FWHMpos_vec = cell2mat(FWHMpos_cell);
    
    figure; ax = axes; hold on;
    colormap(flipud(My_colorMap));
    scatter(Cf_mat(Logical_tot),FWHMall_vec(Logical_tot),[],groupCatagrs(Logical_tot),'s','filled','DisplayName','all slopes');
    scatter(Cf_mat(Logical_tot),FWHMneg_vec(Logical_tot),[],groupCatagrs(Logical_tot),'s','filled','DisplayName','neg slopes');
    scatter(Cf_mat(Logical_tot),FWHMpos_vec(Logical_tot),[],groupCatagrs(Logical_tot),'s','filled','DisplayName','pos slopes');
    title('FWHM vs. C_f');
    ylabel('FWHM [m]');
    xlabel('C_f [m/s]');
end

%% plot kink location
if 0
    rise1XAllcell = arrayfun(@(S) S.rise1XAll,myStruct,'UniformOutput',0);
    rise1XAllcell = rise1XAllcell(:);
    rise1XAll_vec = cell2mat(rise1XAllcell);
    rise1XNegcell = arrayfun(@(S) S.rise1XNeg,myStruct,'UniformOutput',0);
    rise1XNegcell = rise1XNegcell(:);
    rise1XNeg_vec = cell2mat(rise1XNegcell);
    rise1XPoscell = arrayfun(@(S) S.rise1XPos,myStruct,'UniformOutput',0);
    rise1XPoscell = rise1XPoscell(:);
    rise1XPos_vec = cell2mat(rise1XPoscell);
    
    figure; ax = axes; hold on;
    colormap(flipud(My_colorMap));
    scatter(Cf_mat(Logical_tot),rise1XAll_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    scatter(Cf_mat(Logical_tot),rise1XNeg_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','neg slopes');
    scatter(Cf_mat(Logical_tot),rise1XPos_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','pos slopes');
    title('1^{st} kink vs. C_f');
    ylabel('1^{st} kink loc [m]');
    xlabel('C_f [m/s]');
end
%-- second kink
if 0
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
    scatter(Cf_mat(Logical_tot),rise2XAll_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    scatter(Cf_mat(Logical_tot),rise2XNeg_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','neg slopes');
    scatter(Cf_mat(Logical_tot),rise2XPos_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','pos slopes');
    title('2^{nd} kink vs. C_f');
    ylabel('2^{nd} kink loc [m]');
    xlabel('C_f [m/s]');
end
%-- diff between 1 and 2
if 0
    figure; ax = axes; hold on;
    colormap(flipud(My_colorMap));
    scatter(Cf_mat(Logical_tot),rise1XAll_vec(Logical_tot)-rise2XAll_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    scatter(Cf_mat(Logical_tot),rise1XNeg_vec(Logical_tot)-rise2XNeg_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','neg slopes');
    scatter(Cf_mat(Logical_tot),rise1XPos_vec(Logical_tot)-rise2XPos_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','pos slopes');
    
    %     plot(Cf_mat(Logical_tot),rise1XAll_vec(Logical_tot)-rise2XAll_vec(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),rise1XNeg_vec(Logical_tot)-rise2XNeg_vec(Logical_tot),'s','DisplayName',['neg slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),rise1XPos_vec(Logical_tot)-rise2XPos_vec(Logical_tot),'s','DisplayName',['pos slopes ',myStruct(1).Name(1:10)]);
    title('\Deltakink vs. C_f');
    ylabel('\Deltakink [m]');
    xlabel('C_f [m/s]');
    
end

%% v vs. Cf
if 0
    %--- vres
    colormap(flipud(My_colorMap));
    figure; ax = axes; hold on;
    
    scatter(Cf_mat(Logical_tot),V_res_all(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','spatial, all slopes');
    scatter(Cf_mat(Logical_tot),V_res_neg(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','spatial, neg slopes');
    scatter(Cf_mat(Logical_tot),V_res_pos(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','spatial, pos slopes');
    
    %     plot(Cf_mat(Logical_tot),V_res_all(Logical_tot),'*','DisplayName','spatial, all slopes');
    %     plot(Cf_mat(Logical_tot),V_res_neg(Logical_tot),'*','DisplayName','spatial, neg slopes')
    %     plot(Cf_mat(Logical_tot),V_res_pos(Logical_tot),'*','DisplayName','spatial, pos slopes')
    
    xlabel('C_f [m/s]');
    ylabel('v_{res} [m/s]');
    %     ax.XAxis.Scale='log';
    title(['v_{res} vs. C_f, res location=',num2str(VresLoc)]);
    
    
    
    %--- vpeak
    colormap(flipud(My_colorMap));
    figure; ax = axes; hold on;
    
    scatter(Cf_mat(Logical_tot),Vpeak_All(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','spatial, all slopes');
    scatter(Cf_mat(Logical_tot),Vpeak_Neg(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','spatial, neg slopes');
    scatter(Cf_mat(Logical_tot),Vpeak_Pos(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','spatial, pos slopes');
    
    
    %     plot(Cf_mat(Logical_tot),Vpeak_All(Logical_tot),'*','DisplayName','spatial, all slopes');
    %     plot(Cf_mat(Logical_tot),Vpeak_Neg(Logical_tot),'*','DisplayName','spatial, neg slopes')
    %     plot(Cf_mat(Logical_tot),Vpeak_Pos(Logical_tot),'*','DisplayName','spatial, pos slopes')
    
    xlabel('C_f [m/s]');
    ylabel('v_{peak} [m/s]');
    %     ax.XAxis.Scale='log';
    title(['v_{peak} vs. C_f, res location=',num2str(VresLoc)]);
end

%% extract tau_yield
if 0
    yieldSxy1stAll_cell = arrayfun(@(S) S.yieldSxy1stAll,myStruct,'UniformOutput',0);
    yieldSxy1stAll_cell = yieldSxy1stAll_cell(:);
    yieldSxy1stAll = cell2mat(yieldSxy1stAll_cell);
    yieldSxy1stNeg_cell = arrayfun(@(S) S.yieldSxy1stNeg,myStruct,'UniformOutput',0);
    yieldSxy1stNeg_cell = yieldSxy1stNeg_cell(:);
    yieldSxy1stNeg = cell2mat(yieldSxy1stNeg_cell);
    yieldSxy1stPos_cell = arrayfun(@(S) S.yieldSxy1stPos,myStruct,'UniformOutput',0);
    yieldSxy1stPos_cell = yieldSxy1stPos_cell(:);
    yieldSxy1stPos = cell2mat(yieldSxy1stPos_cell);
    
    
    figure; ax = axes; hold on;
    scatter(Cf_mat(Logical_tot),yieldSxy1stAll(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','spatial, all slopes');
    scatter(Cf_mat(Logical_tot),yieldSxy1stNeg(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','spatial, neg slopes');
    scatter(Cf_mat(Logical_tot),yieldSxy1stPos(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','spatial, pos slopes');
    
    %     plot(Cf_mat(Logical_tot),yieldSxy1stAll(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),yieldSxy1stNeg(Logical_tot),'s','DisplayName',['neg slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),yieldSxy1stPos(Logical_tot),'s','DisplayName',['pos slopes ',myStruct(1).Name(1:10)]);
    title('\tau_{xy} @1^{st} kink vs. C_f');
    ylabel('\tau_{yield} [Pa]');
    xlabel('C_f [m/s]');
    
    
end

%% tau yeild scnd
if 0
    yieldSxy2ndAll_cell = arrayfun(@(S) S.yieldSxy2ndAll,myStruct,'UniformOutput',0);
    yieldSxy2ndAll_cell = yieldSxy2ndAll_cell(:);
    yieldSxy2ndAll = cell2mat(yieldSxy2ndAll_cell);
    yieldSxy2ndNeg_cell = arrayfun(@(S) S.yieldSxy2ndNeg,myStruct,'UniformOutput',0);
    yieldSxy2ndNeg_cell = yieldSxy2ndNeg_cell(:);
    yieldSxy2ndNeg = cell2mat(yieldSxy2ndNeg_cell);
    yieldSxy2ndPos_cell = arrayfun(@(S) S.yieldSxy2ndPos,myStruct,'UniformOutput',0);
    yieldSxy2ndPos_cell = yieldSxy2ndPos_cell(:);
    yieldSxy2ndPos = cell2mat(yieldSxy2ndPos_cell);
    
    
    figure; ax = axes; hold on;
    scatter(Cf_mat(Logical_tot),yieldSxy2ndAll(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','spatial, all slopes');
    scatter(Cf_mat(Logical_tot),yieldSxy2ndNeg(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','spatial, neg slopes');
    scatter(Cf_mat(Logical_tot),yieldSxy2ndPos(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','spatial, pos slopes');
    
    %     plot(Cf_mat(Logical_tot),yieldSxy2ndAll(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),yieldSxy2ndNeg(Logical_tot),'s','DisplayName',['neg slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),yieldSxy2ndPos(Logical_tot),'s','DisplayName',['pos slopes ',myStruct(1).Name(1:10)]);
    title('\tau_{xy} @2^{nd} kink vs. C_f');
    ylabel('\tau_{yield} [Pa]');
    xlabel('C_f [m/s]');
end

%% extract residual velocity
if 0
    [~,locLogicV] = min(abs(myStruct(1).residual_locs - VresLoc));
    
    tauResiduals_mat_cell = arrayfun(@(S) S.tauResiduals_mat(locLogicV,:),myStruct,'UniformOutput',0);
    tauResiduals = cell2mat(tauResiduals_mat_cell);
    
    figure; ax = axes; hold on;
    scatter(Cf_mat(Logical_tot),tauResiduals(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o');
    %     plot(Cf_mat(Logical_tot),tauResiduals(Logical_tot),'o','DisplayName',myStruct(1).Name(1:10));
    title('\tau_{res} vs. C_f');
    ylabel('\tau_{res} [MPa]');
    xlabel('C_f [m/s]');
end

%% peak velocity location
if 0
    VpLoc_All_cell = arrayfun(@(S) S.peakLoc_vecAll,myStruct,'UniformOutput',0);
    VpLoc_All_cell = VpLoc_All_cell(:);
    VpLoc_All = cell2mat(VpLoc_All_cell);
    VpLoc_Pos_cell = arrayfun(@(S) S.peakLoc_vecPos,myStruct,'UniformOutput',0);
    VpLoc_Pos_cell = VpLoc_Pos_cell(:);
    VpLoc_Pos = cell2mat(VpLoc_Pos_cell);
    VpLoc_Neg_cell = arrayfun(@(S) S.peakLoc_vecNeg,myStruct,'UniformOutput',0);
    VpLoc_Neg_cell = VpLoc_Neg_cell(:);
    VpLoc_Neg = cell2mat(VpLoc_Neg_cell);
    
    
    figure; ax = axes; hold on;
    scatter(Cf_mat(Logical_tot),VpLoc_All(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','spatial, all slopes');
    scatter(Cf_mat(Logical_tot),VpLoc_Neg(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','spatial, neg slopes');
    scatter(Cf_mat(Logical_tot),VpLoc_Pos(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','spatial, pos slopes');
    %     plot(Cf_mat(Logical_tot),VpLoc_All(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),VpLoc_Neg(Logical_tot),'s','DisplayName',['neg slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),VpLoc_Pos(Logical_tot),'s','DisplayName',['pos slopes ',myStruct(1).Name(1:10)]);
    title('v_{peak,Loc} vs. C_f');
    ylabel('v_{peak,Loc} [m]');
    xlabel('C_f [m/s]');
    
end
%% Cf*VpeakLoc/Vpeak
if 0
    VpLoc_All_cell = arrayfun(@(S) S.peakLoc_vecAll,myStruct,'UniformOutput',0);
    VpLoc_All_cell = VpLoc_All_cell(:);
    VpLoc_All = cell2mat(VpLoc_All_cell);
    VpLoc_Pos_cell = arrayfun(@(S) S.peakLoc_vecPos,myStruct,'UniformOutput',0);
    VpLoc_Pos_cell = VpLoc_Pos_cell(:);
    VpLoc_Pos = cell2mat(VpLoc_Pos_cell);
    VpLoc_Neg_cell = arrayfun(@(S) S.peakLoc_vecNeg,myStruct,'UniformOutput',0);
    VpLoc_Neg_cell = VpLoc_Neg_cell(:);
    VpLoc_Neg = cell2mat(VpLoc_Neg_cell);
    
    
    figure; ax = axes; hold on;
    scatter(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*VpLoc_All(Logical_tot)./Vpeak_All(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    scatter(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*VpLoc_Neg(Logical_tot)./Vpeak_Neg(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','neg slopes');
    scatter(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*VpLoc_Pos(Logical_tot)./Vpeak_Pos(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','pos slopes');
    %     plot(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*VpLoc_All(Logical_tot)./Vpeak_All(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*VpLoc_Neg(Logical_tot)./Vpeak_Neg(Logical_tot),'s','DisplayName',['neg slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*VpLoc_Pos(Logical_tot)./Vpeak_Pos(Logical_tot),'s','DisplayName',['pos slopes ',myStruct(1).Name(1:10)]);
    title('maybe X_c vs. C_f');
    ylabel('C_f*v_{pkLoc}/v_{peak} [m]');
    xlabel('C_f [m/s]');
    
end

%% Xc from kinematics
if 0
    dc = 1.5e-6;
    
    
    figure; ax = axes; hold on;
    scatter(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*dc./Vpeak_All(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    scatter(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*dc./Vpeak_Neg(Logical_tot),[],groupCatagrs(Logical_tot),'filled','d','DisplayName','neg slopes');
    scatter(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*dc./Vpeak_Pos(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','pos slopes');
    
    %     plot(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*dc./Vpeak_All(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*dc./Vpeak_Neg(Logical_tot),'s','DisplayName',['neg slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_mat(Logical_tot),Cf_mat(Logical_tot).*dc./Vpeak_Pos(Logical_tot),'s','DisplayName',['pos slopes ',myStruct(1).Name(1:10)]);
    title('estimate X_c vs. C_f');
    ylabel('C_f*d_{c}/v_{peak} [m]');
    xlabel('C_f [m/s]');
    
end

%% get LEFM distance when v(LEFM)=v(measure)
if 0
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    G_mat = cell2mat(G_Cell);
    
    [~,~, Cr]=CrackSolutionMaterialProperties;
    
    G_filter = G_mat(Logical_tot);
    Cf_filter = Cf_mat(Logical_tot);
    Vpeak_All_filter = Vpeak_All(Logical_tot);
    Vpeak_Neg_filter = Vpeak_Neg(Logical_tot);
    Vpeak_Pos_filter = Vpeak_Pos(Logical_tot);
    N = length(Cf_filter);
    x = -0.05:1e-7:0.005;
    XintersectCell_All = {};
    for i=1:N
        sol = CrackSolutionForh_GammaChange(Cf_filter(i)/Cr , -0.5 , 1e-9 , G_filter(i), x);
        [X0All,~] = intersections(sol.x,sol.vx,[-1 1]*100,[1 1]*Vpeak_All_filter(i));
        [X0Neg,~] = intersections(sol.x,sol.vx,[-1 1]*100,[1 1]*Vpeak_Neg_filter(i));
        [X0Pos,~] = intersections(sol.x,sol.vx,[-1 1]*100,[1 1]*Vpeak_Pos_filter(i));
        XintersectCell_All{i} = X0All;
        XintersectCell_Neg{i} = X0Neg;
        XintersectCell_Pos{i} = X0Pos;
    end
    L = cellfun(@length,XintersectCell_All);
    XintersectCell_filterAll =XintersectCell_All(L~=0);
    G_filter2 = G_filter(L~=0);
    Cf_filter2 = Cf_filter(L~=0);
    Xinter_all = cellfun(@(A) A(1),XintersectCell_filterAll);
    
    L = cellfun(@length,XintersectCell_Neg);
    XintersectCell_filterNeg =XintersectCell_Neg(L~=0);
    G_filter2 = G_filter(L~=0);
    Cf_filter2 = Cf_filter(L~=0);
    Xinter_Neg = cellfun(@(A) A(1),XintersectCell_filterNeg);
    
    L = cellfun(@length,XintersectCell_Pos);
    XintersectCell_filterPos =XintersectCell_Pos(L~=0);
    G_filter2 = G_filter(L~=0);
    Cf_filter2 = Cf_filter(L~=0);
    Xinter_Pos = cellfun(@(A) A(1),XintersectCell_filterPos);
    groupCatagrs_filter = groupCatagrs(Logical_tot);
    
    figure; hold on;
    scatter(Cf_filter2,Xinter_all,[],groupCatagrs_filter(L~=0),'filled','o','DisplayName','all slopes');
    scatter(Cf_filter2,Xinter_Neg,[],groupCatagrs_filter(L~=0),'filled','d','DisplayName','neg slopes');
    scatter(Cf_filter2,Xinter_Pos,[],groupCatagrs_filter(L~=0),'filled','s','DisplayName','pos slopes');
    %     plot(Cf_filter2,Xinter_all,'*','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_filter2,Xinter_Neg,'*','DisplayName',['neg slopes ',myStruct(1).Name(1:10)]);
    %     plot(Cf_filter2,Xinter_Pos,'*','DisplayName',['pos slopes ',myStruct(1).Name(1:10)]);
end
%% errors Vs. Vpeak
if 0
    varVpeak_All_cell = arrayfun(@(S) S.VarVpeak_vecAll,myStruct,'UniformOutput',0);
    varVpeak_All_cell = varVpeak_All_cell(:);
    varVpeak_All = cell2mat(varVpeak_All_cell);
    
    
    figure;
    scatter(Cf_mat(Logical_tot),varVpeak_All(Logical_tot)./Vpeak_All(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o');
    %     plot(Vpeak_All(Logical_tot),varVpeak_All(Logical_tot)./Vpeak_All(Logical_tot),'o');
    
end

%% tau_p based on Cf and Gamma
if 0
    L = 3e-3*10;
    model = 'exp';
    
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    G_mat = cell2mat(G_Cell);
    
    G_filter = G_mat(Logical_tot);
    Cf_filter = Cf_mat(Logical_tot);
    
    N = length(Cf_filter);
    tau_p_vec = nan(N,1);
    for i=1:N
        tau_p_vec(i) = CohesiveSolutionGetTau_p(Cf_filter(i)./Cr,G_filter(i),L,model);
    end
    
    figure;
    scatter(Cf_filter,tau_p_vec,[],groupCatagrs(Logical_tot),'filled','o');
    %     plot(Cf_filter,tau_p_vec,'o');
end

%% diff(x_non_lin,x_peak)
if 0
    rise2XAllcell = arrayfun(@(S) S.rise2XAll,myStruct,'UniformOutput',0);
    rise2XAllcell = rise2XAllcell(:);
    rise2XAll_vec = cell2mat(rise2XAllcell);
    
    VpLoc_All_cell = arrayfun(@(S) S.peakLoc_vecAll,myStruct,'UniformOutput',0);
    VpLoc_All_cell = VpLoc_All_cell(:);
    VpLoc_All = cell2mat(VpLoc_All_cell);
    
    figure; ax = axes; hold on;
    colormap(flipud(My_colorMap));
    %     plot(Cf_mat(Logical_tot),rise2XAll_vec(Logical_tot)-VpLoc_All(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    scatter(Cf_mat(Logical_tot),rise2XAll_vec(Logical_tot)-VpLoc_All(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    title('x_{non-lin}-x_{peak} vs. C_f');
    ylabel('x_{non-lin}-x_{peak} [m]');
    xlabel('C_f [m/s]');
end

%% Arrhenius plots
if 0
    TypLength = 3e-3;
    TypicalTimeVec = TypLength./Cf_mat(Logical_tot);
    figure;
    A = axes;
    T = slmeval(Cf_mat(Logical_tot),slmTvsCf);
    semilogy(1./abs(T),TypicalTimeVec.*yieldSxy2ndAll(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    title('sort of Arrhenius plot');
    ylabel('$\tau_y\cdot\frac{0.003}{C_f} [Pa\cdot s]  (\eta??)$','Interpreter','latex','FontWeight','bold');
    xlabel('1/C_f [m/s]');
end

%% Estimations using linear slip weakening
if 0
    %--- estimate d_c using tau_nl as tau_p
    yieldSxy2ndAll_cell = arrayfun(@(S) S.yieldSxy2ndAll,myStruct,'UniformOutput',0);
    yieldSxy2ndAll_cell = yieldSxy2ndAll_cell(:);
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    
    notEmpty = ~cellfun(@isempty,G_Cell);
    
    yieldSxy2ndAll = cell2mat(yieldSxy2ndAll_cell(notEmpty))/0.02;
    G_mat = cell2mat(G_Cell(notEmpty));
    
    
    d_c = 2*G_mat./yieldSxy2ndAll;
    
    figure;
    scatter(Cf_mat(Logical_tot),d_c(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    %     plot(Cf_mat(Logical_tot),d_c(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    title('d_c calculated as linear slip weak using \tau_{non lin} as tau_p');
    ylabel('d_c [m]');
    xlabel('C_f [m/s]');
    
    
    Xc = 9*pi*G_mat*E./((1-nu^2)*32*yieldSxy2ndAll.^2);
    
    figure;
    scatter(Cf_mat(Logical_tot),Xc(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    %     plot(Cf_mat(Logical_tot),Xc(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    title('X_c calculated as linear slip weak using \tau_{non lin} as tau_p');
    ylabel('X_c [m]');
    xlabel('C_f [m/s]');
    
end

%% estimate d_c by kinematics
if 0
    VpLoc_All_cell = arrayfun(@(S) S.peakLoc_vecAll,myStruct,'UniformOutput',0);
    VpLoc_All_cell = VpLoc_All_cell(:);
    VpLoc_All = cell2mat(VpLoc_All_cell);
    figure;
    scatter(Cf_mat(Logical_tot),Vpeak_All(Logical_tot).*VpLoc_All(Logical_tot)./Cf_mat(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    %     plot(Cf_mat(Logical_tot),Vpeak_All(Logical_tot).*VpLoc_All(Logical_tot)./Cf_mat(Logical_tot),'s','DisplayName',['all slopes ',myStruct(1).Name(1:10)]);
    
end

%% delta_i distributions
if 0
    shifts_vec_cell = arrayfun(@(S) S.shifts_vec,myStruct,'UniformOutput',0);
    shifts_vec_cell2 = vertcat(shifts_vec_cell{:});
    shifts_vec_cellSelected = shifts_vec_cell2(Logical_tot);
    
    shifts_vec_cellSelected2 = vertcat(shifts_vec_cellSelected{:});
    L = cellfun(@length, shifts_vec_cellSelected);
    delta_i_mat = nan(max(L),length(L));
    for i=1:length(L)
        delta_i_mat(1:L(i),i) = shifts_vec_cellSelected{i};
    end
    Cf_filter = Cf_mat(Logical_tot);
    
    
    V = std(delta_i_mat,1,'omitnan');
    figure;
    scatter(Cf_filter,V,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    %     plot(Cf_filter,V,'.');
    %     figure;
    %     boxplot(delta_i_mat,'Positions',Cf_filter,'whisker',inf);
    
    
    %--- estimate tau_p by using delta and linear model:
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    G_mat = cell2mat(G_Cell);
    G_filter = G_mat(Logical_tot);
    
    tau_p_triangleEstimate = 2.*G_filter(:)./mean(delta_i_mat(:),1,'omitnan');
    figure;
    scatter(Cf_filter,tau_p_triangleEstimate,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    %     plot(Cf_filter,tau_p_triangleEstimate,'.');
end

%% tau_p predictor with specific Gamma
if 0
    Model = 'exp';
    
    yieldSxy2ndAll_cell = arrayfun(@(S) S.yieldSxy2ndAll,myStruct,'UniformOutput',0);
    yieldSxy2ndAll_cell = yieldSxy2ndAll_cell(:);
    yieldSxy2ndAll = cell2mat(yieldSxy2ndAll_cell);
    
    %-- gamma vector
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    G_mat = cell2mat(G_Cell);
    G_filter = G_mat(Logical_tot);
    %-- cohesive length value:
    L = 3.2e-3;
    %-- Cf:
    Cf_filter = Cf_mat(Logical_tot);
    %-- peak vel measurement:
    Vpeak_All_filter = Vpeak_All(Logical_tot);
    
    tau_p_vec = nan(size(G_filter));
    peakVel_vec = nan(size(G_filter));
    for i=1:length(G_filter)
        if Cf_filter(i)>=Cr
            thisCf = Cr-1;
        else
            thisCf = Cf_filter(i);
        end
        [X_peak, peakVel_vec(i)] = cohesive_findPeakOfVelByIterations(thisCf/Cr,1e-9,G_filter(i),'L',L*10,Model);
        tau_p_vec(i) = CohesiveSolutionGetTau_p(thisCf/Cr,G_filter(i),L*10,Model);
    end
    
    Ratio = tau_p_vec./peakVel_vec;
    tau_p_semiPredictor = Vpeak_All_filter.*Ratio;
    
    
    slmRatVsCf = cohesive_tau_p_predictorByVx(Model);
    
    Ratio_4Cf = ppval(Cf_filter,slmRatVsCf);
    
    tau_p_byModel = Ratio_4Cf.*Vpeak_All_filter;
    %-- constructing a model:
    TT = log(tau_p_byModel);
    TT_logic = ~isinf(TT);
    log_tau_p_byModel = TT(TT_logic);
    Cf_filter2 = Cf_filter(TT_logic);
    FitResult = fit_piecewiseLinear(Cf_filter2,log_tau_p_byModel,[14.5,1,-1/100,1000],-inf*[1 1 1 1],inf*[1 1 1 1]);
    plusfun = @(x) max(x,0);
    model = @(P,x) P(1) + P(2)*plusfun(P(4)-x) + P(3)*plusfun(x-P(4));
    
    figure; ax = axes; hold on;
    plot(Cf_filter,tau_p_byModel,'ro','DisplayName','simple model (all slopes)');
    semilogy(1:1255,exp(model(FitResult,1:1255)),'b','DisplayName','fit 4 \tau_p by model');    % fitted model
    plot(Cf_mat(Logical_tot),yieldSxy2ndAll(Logical_tot),'go','DisplayName','\tau_{nl} (all slopes)');
    semilogy(FitResult(4),exp(model(FitResult,FitResult(4))),'mo','MarkerSize',12,'DisplayName','kink loc.');  % mark the kink
    text(FitResult(4),exp(model(FitResult,FitResult(4))),['kink loc. ',num2str(round(FitResult(4)))],'HorizontalAlignment','left'); % txt foe the kink marking
    plot(Cf_filter,tau_p_vec,'ko','DisplayName','\tau_p only by C_f');
    plot(Cf_filter,tau_p_semiPredictor,'co','DisplayName','\tau_p by ration w/ v_{peak}');
    
    ax.YAxis.Scale = 'log';
    xlabel('C_f [m/s]');
    ylabel('\tau_p by model');
    title('comparison \tau_p by cohesive model to tau_{tilde}');
    legend('show');
end


%% tau_p predictor - general considiotns (gamma, cohesive length...)
if 0
    Model = 'exp';
    
    yieldSxy2ndAll_cell = arrayfun(@(S) S.yieldSxy2ndAll,myStruct,'UniformOutput',0);
    yieldSxy2ndAll_cell = yieldSxy2ndAll_cell(:);
    yieldSxy2ndAll = cell2mat(yieldSxy2ndAll_cell);
    
    slmRatVsCf = cohesive_tau_p_predictorByVx(Model);
    Cf_filter = Cf_mat(Logical_tot);
    Ratio_4Cf = ppval(Cf_filter,slmRatVsCf);
    Vpeak_All_filter = Vpeak_All(Logical_tot);
    
    tau_p_byModel = Ratio_4Cf.*Vpeak_All_filter;
    %-- constructing a model:
    TT = log(tau_p_byModel);
    TT_logic = ~isinf(TT);
    log_tau_p_byModel = TT(TT_logic);
    Cf_filter2 = Cf_filter(TT_logic);
    FitResult = fit_piecewiseLinear(Cf_filter2,log_tau_p_byModel,[14.5,1,-1/100,1000],-inf*[1 1 1 1],inf*[1 1 1 1]);
    plusfun = @(x) max(x,0);
    model = @(P,x) P(1) + P(2)*plusfun(P(4)-x) + P(3)*plusfun(x-P(4));
    
    figure; ax = axes; hold on;
    scatter(Cf_filter,tau_p_byModel,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    semilogy(1:1255,exp(model(FitResult,1:1255)),'DisplayName','fit 4 \tau_p by model');    % fitted model
    scatter(Cf_mat(Logical_tot),yieldSxy2ndAll(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    semilogy(FitResult(4),exp(model(FitResult,FitResult(4))),'o','MarkerSize',12);  % mark the kink
    text(FitResult(4),exp(model(FitResult,FitResult(4))),['kink loc. ',num2str(round(FitResult(4)))],'HorizontalAlignment','left'); % txt foe the kink marking
    
    %     semilogy(Cf_filter,tau_p_byModel,'o','DisplayName','\tau_p by cohesive model'); % data
    %     semilogy(1:1255,exp(model(FitResult,1:1255)),'DisplayName','fit 4 \tau_p by model');    % fitted model
    %     semilogy(Cf_mat(Logical_tot),yieldSxy2ndAll(Logical_tot),'s','DisplayName','\tau_{tilde,nonLin}');  % tau tilde (non-lin) from measurements
    %     semilogy(FitResult(4),exp(model(FitResult,FitResult(4))),'o','MarkerSize',12);  % mark the kink
    
    ax.YAxis.Scale = 'log';
    xlabel('C_f [m/s]');
    ylabel('\tau_p by model');
    title('comparison \tau_p by cohesive model to tau_{tilde}');
    legend('show');
end
%% cohesive length estimation
if 0
    %-- this gives an ugly result:
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    G_mat = cell2mat(G_Cell);
    G_filter = G_mat(Logical_tot);
    L_estimate = CohesiveSolutionGetL(Cf_filter/Cr,G_filter,tau_p_byModel)/10;    % deviding by 10 is for and exponential model
    %     L_estimate = CohesiveSolutionGetL(Cf_filter/Cr,G_filter,tau_p_byModel)/3;    % deviding by 3 is for the my sigmoid model
    %-- estimation ofr every Cf, fitted tau_p, const Gamma
    L_estimate_model = CohesiveSolutionGetL((1:1255)/Cr,3.5,exp(model(FitResult,1:1255)))/10;     % deviding by 10 is for and exponential model
    %     L_estimate_model = CohesiveSolutionGetL((1:1255)/Cr,3.5,exp(model(FitResult,1:1255)))/3;     % deviding by 3 is for the my sigmoid model
    %--- comparison to X_pk-LEFM
    Vpeak_All_filter = Vpeak_All(Logical_tot);
    N = length(Cf_filter);
    x = -0.05:1e-7:0.005;
    XintersectCell_All = {};
    for i=1:N
        sol = CrackSolutionForh_GammaChange(Cf_filter(i)/Cr , -0.5 , 1e-9 , G_filter(i), x);
        [X0All,~] = intersections(sol.x,sol.vx,[-1 1]*100,[1 1]*Vpeak_All_filter(i));
        XintersectCell_All{i} = X0All;
    end
    L = cellfun(@length,XintersectCell_All);
    XintersectCell_filterAll =XintersectCell_All(L~=0);
    G_filter2 = G_filter(L~=0);
    Cf_filter2 = Cf_filter(L~=0);
    Xinter_all = cellfun(@(A) A(1),XintersectCell_filterAll);
    groupCatagrs_filter = groupCatagrs(Logical_tot);
    
    
    %--- add deviation from LEFM on left:
    maxXi_vel_vec_cell = arrayfun(@(S) S.maxXi_vel_vec_allSlps,myStruct,'UniformOutput',0);
    maxXi_vel_vec_cell = maxXi_vel_vec_cell(:);
    maxXi_vel_vec = cell2mat(maxXi_vel_vec_cell);
    
    maxXi_loc_vec_cell = arrayfun(@(S) S.maxXi_loc_vec_allSlps,myStruct,'UniformOutput',0);
    maxXi_loc_vec_cell = maxXi_loc_vec_cell(:);
    maxXi_loc_vec = cell2mat(maxXi_loc_vec_cell);
    
    Cf_filter = Cf_mat(Logical_tot);
    
    
    
    figure; ax = axes; hold on;
    scatter(Cf_filter,L_estimate,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','cohesive length by \tau_p estimation');
    plot(1:1255,L_estimate_model,'DisplayName','cohesive length per Cf, fitted tau_p, const Gamma');
    scatter(Cf_filter2,abs(Xinter_all),[],groupCatagrs_filter(L~=0),'d','filled','DisplayName','X_{pk-LEFM}');
    scatter(Cf_filter,abs(maxXi_vel_vec(Logical_tot)),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','deviation in velocity');
    xlabel('C_f [m/s]');
    ylabel('L cohesive [m]');
    legend('show');
    
    
    
end

%% total length of deviation
if 0
    rise2XAllcell = arrayfun(@(S) S.rise2XAll,myStruct,'UniformOutput',0);
    rise2XAllcell = rise2XAllcell(:);s
    rise2XAll_vec = cell2mat(rise2XAllcell);
    
    totLength = abs(maxXi_vel_vec(Logical_tot))+abs(rise2XAll_vec(Logical_tot));
    
    figure; ax = axes; hold on;
    colormap(flipud(My_colorMap));
    scatter(Cf_mat(Logical_tot),totLength,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    title('2^{nd} kink vs. C_f');
    ylabel('2^{nd} kink loc [m]');
    xlabel('C_f [m/s]');
end

%% negative x deviation from LEFM
if 0
    maxXi_vel_vec_cell = arrayfun(@(S) S.maxXi_vel_vec,myStruct,'UniformOutput',0);
    maxXi_vel_vec_cell = maxXi_vel_vec_cell(:);
    maxXi_vel_vec = cell2mat(maxXi_vel_vec_cell);
    
    maxXi_loc_vec_cell = arrayfun(@(S) S.maxXi_loc_vec,myStruct,'UniformOutput',0);
    maxXi_loc_vec_cell = maxXi_loc_vec_cell(:);
    maxXi_loc_vec = cell2mat(maxXi_loc_vec_cell);
    
    Cf_filter = Cf_mat(Logical_tot);
    
    figure; hold on;
    scatter(Cf_filter,maxXi_vel_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','deviation in velocity');
    scatter(Cf_filter,maxXi_loc_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','deviation in location');
    
    %     plot(Cf_filter,maxXi_vel_vec(Logical_tot),'o','DisplayName','deviation in velocity');
    %     plot(Cf_filter,maxXi_loc_vec(Logical_tot),'o','DisplayName','deviation in location');
    legend('show');
    
end

%% tau_r/Sigma
if 0
    
    Vpeak_All_filter = Vpeak_All(Logical_tot);
    V_res_all_filter = V_res_all(Logical_tot);
    Cf_filter = Cf_mat(Logical_tot);
    
    tauResiduals_mat_cell = arrayfun(@(S) S.tauResiduals_mat(locLogicV,:),myStruct,'UniformOutput',0);
    tauResiduals = cell2mat(tauResiduals_mat_cell);
    SyyResiduals_mat_cell = arrayfun(@(S) S.SyyResiduals_mat(locLogicV,:),myStruct,'UniformOutput',0);
    SyyResiduals = cell2mat(SyyResiduals_mat_cell);
    SxxResiduals_mat_cell = arrayfun(@(S) S.SyyResiduals_mat(locLogicV,:),myStruct,'UniformOutput',0);
    SxxResiduals = cell2mat(SyyResiduals_mat_cell);
    UxxResiduals_mat_cell = arrayfun(@(S) S.UxxResiduals_mat(locLogicV,:),myStruct,'UniformOutput',0);
    UxxResiduals = cell2mat(UxxResiduals_mat_cell);
    
    v_3_5mm = -Cf_mat(:).*UxxResiduals(:);
    
    figure; hold on;
    scatter(v_3_5mm(Logical_tot),tauResiduals(Logical_tot)./SyyResiduals(Logical_tot),[],groupCatagrs(Logical_tot),'filled','o','DisplayName','deviation in velocity');
    %     scatter(Cf_filter,maxXi_loc_vec(Logical_tot),[],groupCatagrs(Logical_tot),'filled','s','DisplayName','deviation in location');
    
    legend('show');
    
    %--- color by event:
    ii = 2;
    [C,ia,ic] = unique(event_dates);
    V_res_all_filter_date = v_3_5mm(Logical_tot & (ic==ii));
    localFricCoeff = tauResiduals./SyyResiduals;
    localFricCoeff_filter_date = localFricCoeff(Logical_tot & (ic==ii));
    names_date = Names_cell_combined(Logical_tot & (ic==ii));
    [C,ia,exp_groups] = unique(names_date);
    distColors = distinguishable_colors(max(exp_groups));
    groupsColors = distColors(exp_groups,:);
    
    figure; scatter(V_res_all_filter_date,localFricCoeff_filter_date,[],groupsColors,'filled');
    
    CBar = colorbar;
    colormap(distColors);
    CBar.Ticks = movmean(linspace(0,1,max(exp_groups)+1),2,'endpoints','discard');
    [C,ia,exp_groups] = unique(names_date);
    CBar.TickLabels = C;
    names_date = event_details(Logical_tot & (ic==ii));
end

%% local fric coeff vs. res. vel.
if 0
    UxxResiduals_cell = arrayfun(@(S) S.UxxResiduals_mat,myStruct,'UniformOutput',0);
    UxxResiduals_cellAtLoc = cellfun(@(A) squeeze(A(locLogicA,:,:)),UxxResiduals_cell,'UniformOutput',0);
    Cf_cell = arrayfun(@(S) S.cfAtSg,myStruct,'UniformOutput',0);
    Cf_allSg_mat = cell2mat(Cf_cell);
    v_res_cell = cellfun(@(A,B) -A.*B ,UxxResiduals_cellAtLoc ,Cf_cell ,'UniformOutput',0);
    v_res_mat = cell2mat(v_res_cell(:)');
    v_res_mat_filter = v_res_mat(:,Logical_tot);
    
    cf_crit = 0.85*Cr;
    
    SxyResiduals_cell = arrayfun(@(S) S.SxyResiduals_mat,myStruct,'UniformOutput',0);
    SxyResiduals_cellAtLoc = cellfun(@(A) squeeze(A(locLogicA,:,:)),SxyResiduals_cell,'UniformOutput',0);
    SyyResiduals_cell = arrayfun(@(S) S.SyyResiduals_mat,myStruct,'UniformOutput',0);
    SyyResiduals_cellAtLoc = cellfun(@(A) squeeze(A(locLogicA,:,:)),SyyResiduals_cell,'UniformOutput',0);
    localMu_cell = cellfun(@(A,B) A./B ,SxyResiduals_cellAtLoc ,SyyResiduals_cellAtLoc ,'UniformOutput',0);
    localMu_mat = cell2mat(localMu_cell(:)');
    localMu_mat_filter = localMu_mat(:,Logical_tot);
    
    %     figure; plot(v_res_mat_filter(2:14,:),localMu_mat_filter(2:14,:),'.');
    
    %-- groups of experiments
    groupsID = {{'10-2-22','10-20-54'},...
        {'14-28-36','14-48-1','15-4-33'}...
        {'16-50-17','17-18-6','18-11-25'}...
        };
    %     groupsID = {{'10-2-22'},{'10-20-54'},...
    %         {'14-28-36'},{'14-48-1'},{'15-4-33'}...
    %         {'16-50-17'},{'17-18-6'},{'18-11-25'}...
    %         };
    Houres = cellfun(@(S) S((min(strfind(S,' '))+1):end), Names_cell_combined,'UniformOutput',0);
    [C,ia,ic] = unique(event_dates);
    ii = 1;
    figure; hold on;
    %--- color by event:
    for sg_i = [2:4,6:11]
        %-- normalize
        localMu_matNorm = nan(size(localMu_mat));
        for GrpI = 1:length(groupsID)
            fun = @(s)~cellfun('isempty',strfind(Houres,s));
            out = cellfun(fun,groupsID{GrpI},'UniformOutput',0);
            M = ~~sum(cell2mat(out),2);
            Logic_tit_ic_grp = Logical_tot & (ic==ii) & M;
            if nnz(v_res_mat(sg_i,Logic_tit_ic_grp)<=0.1 & v_res_mat(sg_i,Logic_tit_ic_grp)>=0)>0
                meanLowVelVal = mean(localMu_mat(sg_i,Logic_tit_ic_grp&(v_res_mat(sg_i,:)>=0)'&(v_res_mat(sg_i,:)<=0.1)'));
                %                 meanLowVelVal = 1;
                localMu_matNorm(sg_i,Logic_tit_ic_grp) = localMu_mat(sg_i,Logic_tit_ic_grp)./meanLowVelVal;
            else
                localMu_matNorm(sg_i,Logical_tot & (ic==ii)) = nan;
            end
        end
        
        V_res_all_filter_date = v_res_mat(sg_i,Logical_tot & (ic==ii));
        %     localFricCoeff_filter_date = localMu_mat(sg_i,Logical_tot & (ic==ii));
        localFricCoeff_filter_date = localMu_matNorm(sg_i,Logical_tot & (ic==ii));
        relevant_cf_vecs = Cf_allSg_mat(sg_i,Logical_tot & (ic==ii));
        names_date = Names_cell_combined(Logical_tot & (ic==ii));
        
        [C,ia,exp_groups] = unique(names_date);
        distColors = distinguishable_colors(max(exp_groups));
        groupsColors = distColors(exp_groups,:);
        
        lowCfLogic = relevant_cf_vecs<=cf_crit & relevant_cf_vecs>0;
        HighCfLogic = relevant_cf_vecs>cf_crit;
        
        scatter(V_res_all_filter_date(lowCfLogic),localFricCoeff_filter_date(lowCfLogic),[],groupsColors(lowCfLogic(:),:),'o','filled');
        scatter(V_res_all_filter_date(HighCfLogic),localFricCoeff_filter_date(HighCfLogic),[],groupsColors(HighCfLogic(:),:),'p','filled');
    end
    CBar = colorbar;
    colormap(distColors);
    CBar.Ticks = movmean(linspace(0,1,max(exp_groups)+1),2,'endpoints','discard');
    [C,ia,exp_groups] = unique(names_date);
    CBar.TickLabels = C;
    names_date = event_details(Logical_tot & (ic==ii));
    
    
end
%% fracture energy
if 0
    
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    G_mat = cell2mat(G_Cell);
    G_filter = G_mat(Logical_tot);
    
    Cf_filter = Cf_mat(Logical_tot);
    figure; ax = axes; hold on;
%     scatter(Cf_filter,G_filter,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','\Gamma vs. C_f');
    group1 = groupCatagrs(Logical_tot)==1;
    group2 = groupCatagrs(Logical_tot)==2;
    scatter(Cf_filter(group1)/Cr,G_filter(group1),'filled','o','DisplayName','2018-08-29');
    scatter(Cf_filter(group2)/Cr,G_filter(group2),'filled','o','DisplayName','2018-10-10');
end

%% compare residual velocirt from phedi and sg
if 0
    
    groupCatagrs_filter = groupCatagrs(Logical_tot);
    V_res_all_filter = V_res_all(Logical_tot);
    Cf_filter = Cf_mat(Logical_tot);
    
    UxxResiduals_mat_cell = arrayfun(@(S) S.UxxResiduals_mat(locLogicV,8,:),myStruct,'UniformOutput',0);
    UxxResiduals_mat_cellLin = cellfun(@(A) A(:)', UxxResiduals_mat_cell,'UniformOutput',0);
    UxxResiduals = cell2mat(UxxResiduals_mat_cellLin);
    
    v_3_5mm = -Cf_mat(:).*UxxResiduals(:);
    v_3_5mm_filter = v_3_5mm(Logical_tot);
    
    figure; hold on;
    scatter(v_3_5mm_filter(groupCatagrs_filter==1),V_res_all_filter(groupCatagrs_filter==1),[],abs(Cf_filter(groupCatagrs_filter==1)),'o','filled');
    scatter(v_3_5mm_filter(groupCatagrs_filter==2),V_res_all_filter(groupCatagrs_filter==2),[],abs(Cf_filter(groupCatagrs_filter==2)),'s','filled');
    xlabel('v_{res} @ 3.5mm');
    ylabel('v_{res} @ interface');
    hold on;
    axis equal
    axis square
    plot([0 0.3],[0 0.3],'k');
    colorbar;
end

%% compare parameters to estimate FE
if 0
    
    group1_filter = groupCatagrs(Logical_tot)==1;
    group2_filter = groupCatagrs(Logical_tot)==2;
    
    
    Nsg_vec = [6,8,11];
    Cf_filter = Cf_mat(Logical_tot);
    %--- Uxx max
    figure; hold on;
    for Nsg = Nsg_vec
        UxxMaxAmpCell = arrayfun(@(A) squeeze(A.UxxMaxAmpAndLoc(1,Nsg,:)),myStruct,'UniformOutput',0);
        UxxMaxAmp = cell2mat(UxxMaxAmpCell(:));
        UxxMaxAmp_filter = UxxMaxAmp(Logical_tot);
        plot(Cf_filter,UxxMaxAmp_filter,'o','DisplayName',num2str(Nsg));
    end
    title('Uxx vs. C_f')
    minimizingGamma = gamma_findGammaFit_byUijAmp(Cf_filter(group1_filter),UxxMaxAmp_filter(group1_filter));
    [UxxAmp,UyyAmp,UxyAmp,Cf_vec] = CrackStrainAmpVsCf(minimizingGamma);
    hold on; plot(Cf_vec,UxxAmp,'DisplayName',['G=',num2str(minimizingGamma)]);
    
    %     %--- Uyy max
    %     figure; hold on;
    %     for Nsg = Nsg_vec
    %         UyyMaxAmpCell = arrayfun(@(A) squeeze(A.UyyMaxAmpAndLoc(1,Nsg,:)),myStruct,'UniformOutput',0);
    %         UyyMaxAmp = cell2mat(UyyMaxAmpCell(:));
    %         UyyMaxAmp_filter = UyyMaxAmp(Logical_tot);
    %         plot(Cf_filter,UyyMaxAmp_filter,'o','DisplayName',num2str(Nsg));
    %     end
    %     title('Uyy vs. C_f')
    %
    %     %--- Uxy max
    %     figure; hold on;
    %     for Nsg = Nsg_vec
    %         UxyMaxAmpCell = arrayfun(@(A) squeeze(A.UxyMaxAmpAndLoc(1,Nsg,:)),myStruct,'UniformOutput',0);
    %         UxyMaxAmp = cell2mat(UxyMaxAmpCell(:));
    %         UxyMaxAmp_filter = UxyMaxAmp(Logical_tot);
    %         plot(Cf_filter,UxyMaxAmp_filter,'o','DisplayName',num2str(Nsg));
    %     end
    %     title('Uxy vs. C_f')
end

%% get local gamma distribution
if 0
    locGamma_vec_cell = arrayfun(@(S) S.localGvec_vec,myStruct,'UniformOutput',0);
    locGamma_vec_cell2 = vertcat(locGamma_vec_cell{:});
    locGamma_vec_cellSelected = locGamma_vec_cell2(Logical_tot);
    
    shifts_vec_cellSelected2 = vertcat(locGamma_vec_cellSelected{:});
    L = cellfun(@length, locGamma_vec_cellSelected);
    locGamma_i_mat = nan(max(L),length(L));
    for i=1:length(L)
        locGamma_i_mat(1:L(i),i) = locGamma_vec_cellSelected{i};
    end
    Cf_filter = Cf_mat(Logical_tot);
    
    
    V = mean(locGamma_i_mat,1,'omitnan');
    figure;
    scatter(Cf_filter,V,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    %     plot(Cf_filter,V,'.');
    %     figure;
    %     boxplot(delta_i_mat,'Positions',Cf_filter,'whisker',inf);
    
    
    %--- estimate tau_p by using delta and linear model:
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    G_mat = cell2mat(G_Cell);
    G_filter = G_mat(Logical_tot);
    
    tau_p_triangleEstimate = 2.*G_filter(:)./mean(locGamma_i_mat(:),1,'omitnan');
    figure;
    scatter(Cf_filter,tau_p_triangleEstimate,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    
end

%% comapre interface Gamma to Gamma from various Sgs

if 0
    meanLocG_vec = nan(length(myStruct),1);
    for i=1:length(myStruct)
         meanLocG_vec(i) = cellfun(@mean, myStruct(i).localGvec_vec);
    end
    
    
    
    locGamma_vec_cell = arrayfun(@(S) S.localGvec_vec,myStruct,'UniformOutput',0);
    locGamma_vec_cell2 = vertcat(locGamma_vec_cell{:});
    locGamma_vec_cellSelected = locGamma_vec_cell2(Logical_tot);
    
    shifts_vec_cellSelected2 = vertcat(locGamma_vec_cellSelected{:});
    L = cellfun(@length, locGamma_vec_cellSelected);
    locGamma_i_mat = nan(max(L),length(L));
    for i=1:length(L)
        locGamma_i_mat(1:L(i),i) = locGamma_vec_cellSelected{i};
    end
    Cf_filter = Cf_mat(Logical_tot);
    
    
    V = mean(locGamma_i_mat,1,'omitnan');
    figure;
    scatter(Cf_filter,V,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    %     plot(Cf_filter,V,'.');
    %     figure;
    %     boxplot(delta_i_mat,'Positions',Cf_filter,'whisker',inf);
    
    
    %--- estimate tau_p by using delta and linear model:
    G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
    G_Cell = G_Cell(:);
    G_mat = cell2mat(G_Cell);
    G_filter = G_mat(Logical_tot);
    
    tau_p_triangleEstimate = 2.*G_filter(:)./mean(locGamma_i_mat(:),1,'omitnan');
    figure;
    scatter(Cf_filter,tau_p_triangleEstimate,[],groupCatagrs(Logical_tot),'filled','o','DisplayName','all slopes');
    
end
