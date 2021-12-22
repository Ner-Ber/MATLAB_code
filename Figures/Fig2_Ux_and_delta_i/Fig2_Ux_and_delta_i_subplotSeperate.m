%% create figure of system

%% load for sasmples
% load('G:\Frics\2018-10-10\allPhediStructures_17186.mat');     % old
load('G:\Frics\2018-10-10\allPhediStructures_165017.mat');	% new, changed on 2019-06-27
[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
for iEv = 1:length(relevantEvents)
    PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
    PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
end

EvNum = 21;


%% displacement and delta_i
Ux_ax_F = figure;
% Movie_phedi_plotWithSg(PhediStructCell{EvNum},'00',{'locSpace'},[],[],[],Ux_ax_F);
PhediData = PhediStructCell{EvNum}.PhediData;
%-- plot colors
Nphd = size(PhediData.PhediLocation,2);       % total number of phedis
FigColors = MyVaryColor(Nphd);
%-- plot legend
LGND = cellfun(@num2str,num2cell(PhediData.measuredPhedisFromPlot),'UniformOutput',0);
%-- reduce sin
PhediLocReduced = phedi_reduceSinFromLocation(PhediData);
PhediLocReduced = bsxfun(@minus,PhediLocReduced,mean(PhediLocReduced(1:30,:),'omitnan'));

hold on;

%--- set spatial axis
spaceAxis = PhediData.x_mins_x_tip;

%-- set limits for relevant plotting
inBlockLogicals = logical(PhediData.inBlockLogicals);

%--- plot it:
for i = 1:Nphd
    plot(spaceAxis(inBlockLogicals(:,i),i),PhediLocReduced(inBlockLogicals(:,i),i),...
        '-','LineWidth',2,'MarkerSize',4,'Color',FigColors(i,:),'DisplayName',LGND{i});
end
hold off;

% Ux_ax = get(Ux_ax_F,'children');
% for i=1:length(Ux_ax.Children)
%     Ux_ax.Children(i).Marker = '.';
%     Ux_ax.Children(i).MarkerSize = 4;
%     Ux_ax.Children(i).LineWidth = 2;
% end
% Ux_ax.Title.Visible = 'off';
xlabel('x-x_{tip} [m]');
xlim([-0.05 0.05]);
ylabel('u_x [\mum]');
ylim([-1 12]*1e-6);
Ux_ax_F.Children.YAxis.TickValues = [0 5 10]*1e-6;
Ux_ax_F.Children.YAxis.TickLabels = [0 5 10];

% hold on;
% plot(PhediStructCell{EvNum}.solAtInter.x,PhediStructCell{EvNum}.solAtInter.Ux,'k','LineWidth',3);

Ux_ax_F.Name = 'phedi_Ux';
%% phedi histogram
[residualsCell, RMS_vec, stretch_vec] = phedi_UxBestFitVerticalStretch(PhediStructCell{EvNum},[],0);


fixedK = PhediStructCell{EvNum}.solAtSG.K*stretch_vec;
v = PhediStructCell{EvNum}.PhediData.Cf;
k=Cs/Cd;%Broberg p.330
alpha_d=(1-(v./Cd).^2).^0.5;
alpha_s=(1-(v./Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
%----- Following Broberg p.334,336
A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2;  % known as Y_II in broberg p.336
localG = A.*(fixedK.^2)./(mu*4*(1-k^2));

Gloc_hist_F = figure;
histogram(localG);
Gloc_hist = Gloc_hist_F.Children;
hold on;
%-- add G from sg marking
YLim = Gloc_hist.YAxis.Limits;
Gsg = PhediStructCell{EvNum}.solAtSG.Gamma;
plot([1 1]*Gsg,YLim,'LineWidth',1.5);
ylabel('count');
xlabel('\Gamma_{loc} [Jm^{-2}]');

Gloc_hist_F.Name = 'Gloc_hist';

%% plot average phedi
UxMean_ax_F = figure;

AvgPhediStruct = phedi_averagePhedi(PhediStructCell{EvNum});
hold on;
%--- plot errs cloud
XX = [AvgPhediStruct.X_mean;flipud(AvgPhediStruct.X_mean)];
YY = [AvgPhediStruct.Loc_mean_x-sqrt(AvgPhediStruct.Loc_var_x);flipud(AvgPhediStruct.Loc_mean_x+sqrt(AvgPhediStruct.Loc_var_x))];
YY = YY(~isnan(XX));    %-- eliminate Nan
XX = XX(~isnan(XX));
XX = XX(~isnan(YY));
YY = YY(~isnan(YY));
% ErrCloud = fill(XX,YY,rgb('PeachPuff'),'linestyle','none','FaceAlpha',0.5,'DisplayName','phedi Vel std');

%--- Add Err Bars
ErrBarLocs = [0.005 -0.005 -0.02 -0.04];
[~,ErrBarCoor] = min(abs(bsxfun(@minus, AvgPhediStruct.X_mean(:),ErrBarLocs(:)')));
Err = sqrt(AvgPhediStruct.Loc_var_x(ErrBarCoor));
ErrBars = errorbar(AvgPhediStruct.X_mean(ErrBarCoor),AvgPhediStruct.Loc_mean_x(ErrBarCoor),Err,'.','Color',rgb('PeachPuff'),'LineWidth',2);
%-- plot mean phedi
meanCurve = plot(AvgPhediStruct.X_mean,AvgPhediStruct.Loc_mean_x,'-','Color',rgb('Crimson'),'LineWidth',3,'DisplayName','mean phedi');

%-- plot LEFM on interface
LEFMpred = plot(PhediStructCell{EvNum}.solAtInter.x,PhediStructCell{EvNum}.solAtInter.Ux,'k-','LineWidth',2,'DisplayName','LEFM @intrface');


%--- make plot pretty
xlim([-0.05 0.05]);
ylim([-0.02*1e-6 1.3*max(YY(abs(XX)<=0.05))])
xlabel('x-x_{tip} [m]');
ylabel('u_x [m/s]');
hold off;
UxMean_ax_F.Name = 'meanUxPhedi';

%% plot sg sample as inset
sg_sample_F = figure;
Nsg = 8;
hold on;
%--- plot Uxx
plot(PhediStructCell{EvNum}.SgData.x_mins_x_tips(:,Nsg),1e-3*(PhediStructCell{EvNum}.SgData.Uxx(:,Nsg)-mean(PhediStructCell{EvNum}.SgData.Uxx(1:100,Nsg))),...
    '.-','Color',rgb('DimGray'),'LineWidth',1,'MarkerSize',4);
plot(PhediStructCell{EvNum}.solAtSG.x,PhediStructCell{EvNum}.solAtSG.Uxx,...
    'Color',rgb('Chocolate'));
ylim([-2.5 0.5]*1e-4);
xlim([-0.05 0.05]*1e-0);
sg_sample = sg_sample_F.Children;
sg_sample.YAxis.TickValues = [-2 -1 0]*1e-4;
sg_sample.YAxis.TickLabels = [-20 -10 0];
xlabel('x-x_{tip} [m]');
ylabel('\Delta\epsilon_{xx} [10^{-3}]');
xlim([-0.05 0.05]);

sg_sample_F.Name = 'sgSample_inset';