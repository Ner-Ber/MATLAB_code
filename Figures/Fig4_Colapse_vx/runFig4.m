%% create figure of system
[Cd, Cs, Cr, nu, ro, E, mu, ~,~,~,~]=CrackSolutionMaterialProperties;
smthWindow = 9;
maxSrchRegn = 0.03;
Gamma = 3.5;
L = 10*3.2e-3; %3*3.2e-3;
model = 'exp'; %'MySigmoid';
x_vec = -0.03:1e-5:0.01;

%% load for sasmples
% load('G:\Frics\2018-10-10\allPhediStructures_17186.mat');
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
for iEv = 1:length(relevantEvents)
    PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
    PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
end

FontSize = 14;

displayLowEv = [6, 9, 15];
displayHghEv = [12, 16, 17];
displayCmprEv = [9, 16];
displayCmpr2Mdl = [15, 16];

%% create figure and axes
Dx = 1/22; dx = 1/66; Wdth = 3/12; lftSdXtra = 3*dx;
yLevl = 0.2; yHght = 0.7;


F4 = figure;
F4.WindowStyle = 'normal';
F4.Units = 'centimeters';
F4.Position([3,4]) = [32 10];

%-- low Cf axes
vx_lowCfClps = axes;
lowLeftSd = Dx+lftSdXtra;
vx_lowCfClps.Position = [lowLeftSd, yLevl, Wdth, yHght];    % low Cf collapsed

%-- high Cf axes
vx_highCfClps = axes;
hghLftSd = lftSdXtra+Dx+Wdth+dx;
vx_highCfClps.Position = [hghLftSd, yLevl, Wdth, yHght];    % high Cf collapsed
vx_highCfClps.YAxis.TickLabels=[];

%-- CB and compare axes
vx_CmprCf = axes;

CBar = colorbar;
TickVals = [0 0.25 0.5 0.75 1]*Cr/1260;
CBar.Ticks = TickVals;
CBar.TickLabels = cat(2,{0 0.25 0.5 0.75},'C_R');
CBar.Label.String = 'C_f/C_R ';
cbLftSd = lftSdXtra+Dx+3*Wdth+3*dx;
CBar.Position = [cbLftSd, yLevl, 0.5*Dx, yHght];   % colorbar position

cmbnLftSd = lftSdXtra+Dx+2*Wdth+2*dx;
vx_CmprCf.Position = [cmbnLftSd, yLevl, Wdth, yHght];   % comparison between Cfs
vx_CmprCf.YAxis.TickLabels=[];

%-- inset
instWdth = 0.4*Wdth; instHght = 0.33*yHght;
isntdx = 2*dx; isntdy = 2*dx;
Cmpr2Mdl = axes;
Cmpr2Mdl.Position = [cmbnLftSd+isntdx, yLevl+yHght-(isntdy+instHght), instWdth, instHght];    % high Cf collapsed




%% plot normalized curves in all 3 axes
Cf_vec = [];
for iEv = 1:length(relevantEvents)
    Cf_vec(iEv) = PhediStructCell{relevantEvents(iEv)}.PhediData.Cf;
end
relevantEvents(Cf_vec>1260) = [];
Cf_vec(Cf_vec>1260) = [];

relevant_colorscale = linspace(0,Cr,1e3);
relevantColorMap = colormap(MyVaryColor(1e3,flipud(My_colorMap)));
for iEv = 1:length(relevantEvents)
    Name = [PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details,' Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)];
    AvgPhediStruct = phedi_averagePhedi(PhediStructCell{relevantEvents(iEv)});
    x_Phedi = AvgPhediStruct.X_mean;
    vel_phedi = AvgPhediStruct.Vel_mean_x;
    [~,I] = min(abs(Cf_vec(iEv)-relevant_colorscale));
    
    %-- smooth and normalize
    Y = smooth(vel_phedi,smthWindow);
    vel_phedi_norm = Y/max(Y(abs(x_Phedi)<=maxSrchRegn));
    
    %-- select only data before relected waves
%     [~,returnIdx] = min(abs(AvgPhediStruct.T_mean-0.14/Cd));
%     returnLogical = x_Phedi>=AvgPhediStruct.X_mean(returnIdx);
%     [x_Phedi_crp,vel_phedi_norm] = deal(x_Phedi(returnLogical),vel_phedi_norm(returnLogical));
    
    
    axes(vx_lowCfClps); hold on;
    plot(x_Phedi,vel_phedi_norm,...
        '.-','Color',relevantColorMap(I,:),...
        'DisplayName',[PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details, '  Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)]);
    
    axes(vx_highCfClps); hold on;
    plot(x_Phedi,vel_phedi_norm,...
        '.-','Color',relevantColorMap(I,:),...
        'DisplayName',[PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details, '  Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)]);
    
    axes(vx_CmprCf); hold on;
    plot(x_Phedi,vel_phedi_norm,...
        '.-','Color',relevantColorMap(I,:),...
        'DisplayName',[PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details, '  Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)]);
    
    if ismember(relevantEvents(iEv),displayCmpr2Mdl);
        axes(Cmpr2Mdl); hold on;
        plot(x_Phedi,vel_phedi,...
            '.-','Color',relevantColorMap(I,:),...
            'DisplayName',[PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.details, '  Cf=',num2str(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf)]);
        
        %     sol_cohesive = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_vec(iEv)/Cr,1e-9,ppval(GammaFit,Cf_vec(iEv)),'L',L,model,[],x_vec);
        sol_cohesive = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_vec(iEv)/Cr,1e-9,Gamma,'L',L,model,[],x_vec);
        plot(sol_cohesive.x,sol_cohesive.vx,'--','LineWidth',2,'Color',relevantColorMap(I,:));
    end
end

%% select curves for low Cf collapse

flippedEv = rot90(relevantEvents,90);
for i=1:length(vx_lowCfClps.Children)
    if ismember(flippedEv(i),displayLowEv)
        vx_lowCfClps.Children(i).Visible = 'on';
        continue
    end
    vx_lowCfClps.Children(i).Visible = 'off';
end
axes(vx_lowCfClps)
colormap(relevantColorMap);
xlim([-0.055 0.055]*1e-0);
ylim([-0.05 1.3]*1e-0);
xlabel('x-x_{tip} [m]');
ylabel('v_{x, normalized} [m/s]');


%% select curves for high Cf collapse

flippedEv = rot90(relevantEvents,90);
for i=1:length(vx_highCfClps.Children)
    if ismember(flippedEv(i),displayHghEv)
        vx_highCfClps.Children(i).Visible = 'on';
        continue
    end
    vx_highCfClps.Children(i).Visible = 'off';
end

axes(vx_highCfClps)
colormap(relevantColorMap);
xlim([-0.055 0.055]*1e-0);
ylim([-0.05 1.3]*1e-0);
xlabel('x-x_{tip} [m]');


%% high and low Cf comparison

flippedEv = rot90(relevantEvents,90);
for i=1:length(vx_highCfClps.Children)
    if ismember(flippedEv(i),displayCmprEv)
        vx_CmprCf.Children(i).Visible = 'on';
        continue
    end
    vx_CmprCf.Children(i).Visible = 'off';
end

axes(vx_CmprCf)
colormap(relevantColorMap);
xlim([-0.055 0.055]*1e-0);
ylim([-0.05 1.3]*1e-0);
xlabel('x-x_{tip} [m]');


%% compare 2 model inset
axes(Cmpr2Mdl)
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

