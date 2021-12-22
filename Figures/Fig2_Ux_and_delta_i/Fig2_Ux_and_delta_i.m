%% create figure of system

%% load for sasmples
% load('G:\Frics\2018-10-10\allPhediStructures_17186.mat');     % old
% load('G:\Frics\2018-10-10\allPhediStructures_165017.mat');	% new, changed on 2019-06-27
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
for iEv = 1:length(relevantEvents)
    PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
    PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
end

EvNum = 21;
FontSize = 10;

%% create figure and axes
F2 = figure;
F2.Color=[1 1 1];
F2.WindowStyle = 'normal';
F2.Units = 'centimeters';
F2.Position([3,4]) = [18 8];
Ux_ax = axes;
Ux_ax.Position = [0.13, 0.6, 0.35, 0.35];   % top left
delta_i_hist = axes;
delta_i_hist.Position = [0.58, 0.6, 0.35, 0.35];    % top right
UxMean_ax = axes;
UxMean_ax.Position = [0.13, 0.12, 0.8, 0.35];   % bottom
sg_sample = axes;
sg_sample.Position = [0.7, 0.25, 0.2, 0.2];   % inset

%% displacement and delta_i
% Ux_fig = figure;
Movie_phedi_plotWithSg(PhediStructCell{EvNum},'00',{'locSpace'},[],[],[],Ux_ax);
for i=1:length(Ux_ax.Children)
    Ux_ax.Children(i).Marker = '.';
    Ux_ax.Children(i).MarkerSize = 4;
    Ux_ax.Children(i).LineWidth = 2;
end
Ux_ax.Title.Visible = 'off';
ylabel('u_x [\mum]');
ylim([-1 12]*1e-6);
Ux_ax.YAxis.TickValues = [0 5 10]*1e-6;
Ux_ax.YAxis.TickLabels = [0 5 10];

hold on;
% plot(PhediStructCell{EvNum}.solAtInter.x,PhediStructCell{EvNum}.solAtInter.Ux,'k','LineWidth',3);

%--- plot collapsed phedis after delta reduced
% [residualsCell, RMS_vec, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStructCell{EvNum},[],0);
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

slctdPhedi4Exmpl = 15;
di = stretch_vec(slctdPhedi4Exmpl);
%-- create shape of error bar
earL = 0.005;
% delta_y = [di/2 di/2 di/2 -di/2 -di/2 -di/2];
% delta_x = [-earL/2 earL/2 0 0 earL/2 -earL/2];
delta_y = [di/2 -di/2];
delta_x = [0 0];
mark_x = -0.052;
mark_y = 9e-6;




%% phedi histogram
axes(delta_i_hist);
histogram(localG);
hold on;
%-- add G from sg marking
YLim = delta_i_hist.YAxis.Limits;
Gsg = PhediStructCell{EvNum}.solAtSG.Gamma;
plot([1 1]*Gsg,YLim,'LineWidth',1.5);
ylabel('count');
xlabel('\delta_i [\mum]');


%% plot average phedi
axes(UxMean_ax);
Movie_phedi_plotWithSg(PhediStructCell{EvNum},'00',{'avgPhediLoc'},[],[],[],UxMean_ax);
UxMean_ax.Title.Visible = 'off';
UxMean_ax.Children(4).Visible = 'off';
UxMean_ax.Children(5).Visible = 'off';
UxMean_ax.Children(1).Color = 'k';

%% plot sg sample as inset
axes(sg_sample);
Nsg = 8;
hold on;
%--- plot Uxx
plot(PhediStructCell{EvNum}.SgData.x_mins_x_tips(:,Nsg),1e-3*(PhediStructCell{EvNum}.SgData.Uxx(:,Nsg)-mean(PhediStructCell{EvNum}.SgData.Uxx(1:100,Nsg))),...
    '.-','Color',rgb('DimGray'),'LineWidth',1,'MarkerSize',4);
plot(PhediStructCell{EvNum}.solAtSG.x,PhediStructCell{EvNum}.solAtSG.Uxx,...
    'Color',rgb('Chocolate'));
ylim([-2.5 0.5]*1e-4);
xlim([-0.05 0.05]*1e-0);

sg_sample.YAxis.TickValues = [-2 -1 0]*1e-4;
sg_sample.YAxis.TickLabels = [-20 -10 0];
xlabel('x-x_{tip} [m]');
ylabel('\Delta\epsilon_{xx} [10^{-3}]');
xlim([-0.05 0.05]);


%% set all sub figure
set(findall(F2,'-property','FontSize'),'FontSize',FontSize);