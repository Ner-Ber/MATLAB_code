%% load phedi struct
% load('E:\Frics\2018-10-10\allPhediStructures_17186.mat');
[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
for iEv = 1:length(relevantEvents)
    PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
    PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
end

%% load accumulated data
experName = '2018-10-10 17-18-6';
A1 = load('E:\Frics\2018-10-10\DAandREsidualsStruct.mat');
myStruct = A1.myStruct;
%-- take relevant experiment
relevantExper = find(arrayfun(@(S) strcmp(S.Name,experName),myStruct));

Events = myStruct(relevantExper).Events;
X_pk = myStruct(relevantExper).peakLoc_vecAll;
X_nl = myStruct(relevantExper).rise1XAll;

%% set some parameters
optionalEvents = [12, 15, 21];
EvNum = 12;

retrnWave4Cutoff = Cr;
Nsg = 8;
smthWindow = 65;
maxSrchRegn = 0.03;
L = 10*3.2e-3; %3*3.2e-3;
model = 'exp'; %'MySigmoid';
Gamma = 3.5;

shift_vals = linspace(0,4.5e-3, 5);
% shiftColor = colormap(MyVaryColor(length(shift_vals),flipud(My_colorMap)));
% shiftColor = colormap(MyVaryColor(3,flipud(My_colorMap)));
shiftColor = [
    0, 114, 189;
    237, 177, 32;
    53, 155, 22
        ]./255;

X_pk_meas = [-0.00133   -0.00601   -0.00312];


%% plot figures for comparison
for EvNum = optionalEvents
    
    %% compare to cohesive model
    x_vec = -0.06:1e-4:0.02;
    figure('Name',sprintf('Ev%d_phedi',EvNum));
    Cmpr2Mdl = axes; hold on;
    % for iEv = 1:length(Cf_vec)
    AvgPhediStruct = phedi_averagePhedi(PhediStructCell{EvNum});
    yeildSall = yield_findYieldStressStrain(PhediStructCell{EvNum},0.01);
    %     Sall = phedi_FWHMofMean(PhediStructCell{EvNum},0.01);
    %     Sall = phedi_FWHMofMean(PhediStructCell{EvNum},0.002);
    x_Phedi = AvgPhediStruct.X_mean;
    vel_phedi = AvgPhediStruct.Vel_mean_x;
    
    %-- smooth and normalize
    Y = smooth(vel_phedi,smthWindow);
    
    %-- select only data before relected waves
    [returnLogical, returnIdx] = cutSignalsFromReturningWaves(retrnWave4Cutoff,PhediStructCell{EvNum}, AvgPhediStruct);
    %     [x_Phedi,vel_phedi] = deal(x_Phedi(returnLogical),vel_phedi(returnLogical));
    [x_Phedi,vel_phedi] = deal(x_Phedi(returnLogical),Y(returnLogical));
    
    plot(x_Phedi,vel_phedi,...
        '-','Color',[1 1 1]*0.24,...
        'DisplayName',[PhediStructCell{EvNum}.BigPicRotStruct.details,sprintf('  Cf=%0.0f',PhediStructCell{EvNum}.PhediData.Cf)]);
    
    %     sol_cohesive = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_vec(iEv)/Cr,1e-9,ppval(GammaFit,Cf_vec(iEv)),'L',L,model,[],x_vec);
    Cf_ev = PhediStructCell{EvNum}.PhediData.Cf;
    sol_cohesive = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_ev/Cr,1e-9,Gamma,'L',L,model,[],x_vec);
    
    [~,I] = max(sol_cohesive.vx);
    cohesive_pk = sol_cohesive.x(I);
    %     shift_vals = [0, X_pk(Events==EvNum)-cohesive_pk, yeildSall.VelKinkX];
    %     shift_vals = [0, Sall.peakX-cohesive_pk, yeildSall.VelKinkX];
    shift_vals = [0, X_pk_meas(optionalEvents==EvNum)-cohesive_pk, yeildSall.VelKinkX];
    %     shift_vals = [0, X_pk(Events==EvNum)-cohesive_pk, X_nl(Events==EvNum), yeildSall.VelKinkX];
    for l = 1:length(shift_vals)
        plot(sol_cohesive.x+shift_vals(l),sol_cohesive.vx,'--','LineWidth',2, 'Color', shiftColor(l,:), 'DisplayName', num2str(shift_vals(l)));
    end
    
    
    Cmpr2Mdl.FontSize = 10;
    Cmpr2Mdl.YAxis.TickValues = [0 0.2 0.4];
    Cmpr2Mdl.XAxis.TickValues = [-0.05 -0.02 0 0.02];
    xlim([-0.05 0.02]);
    ylim([-0.02 0.55]);
    xlabel('x-x_{tip} [m]');
    ylabel('v_x [m/s]');
    pbaspect([8,5,1])
    legend('show', 'Location','northwest')
    
    
    %% plot sg Uxx fit:
    %     sg_sample_F = figure;
    sg_sample_F = figure('Name',sprintf('Ev%d_sg',EvNum));
    hold on;
    plot(PhediStructCell{EvNum}.SgData.x_mins_x_tips(:,Nsg),1e-3*(PhediStructCell{EvNum}.SgData.Uxx(:,Nsg)-mean(PhediStructCell{EvNum}.SgData.Uxx(1:100,Nsg))),...
        '-','Color',rgb('DimGray'),'LineWidth',1,'MarkerSize',4,...
        'DisplayName',[PhediStructCell{EvNum}.BigPicRotStruct.details, sprintf('  Cf=%0.0f',PhediStructCell{EvNum}.PhediData.Cf)]);
    %     plot(PhediStructCell{EvNum}.solAtSG.x,PhediStructCell{EvNum}.solAtSG.Uxx,...
    %         'Color',rgb('Chocolate'), 'DisplayName', 'LEFM');
    
    sol_cohesive_atSg = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_ev/Cr,3.5e-3,Gamma,'L',L,model,[],x_vec);
    
    for l = 1:length(shift_vals)
        plot(PhediStructCell{EvNum}.solAtSG.x+shift_vals(l),PhediStructCell{EvNum}.solAtSG.Uxx,...
            'Color',shiftColor(l,:)*0.8, 'DisplayName', sprintf('LEFM+%d',shift_vals(l)));
        plot(sol_cohesive_atSg.x+shift_vals(l),sol_cohesive_atSg.Uxx,'--','LineWidth',2, 'Color', shiftColor(l,:), 'DisplayName', num2str(shift_vals(l)));
    end
    
    ylim([-2.5 0.5]*1e-4);
    xlim([-0.035 0.035]);
    sg_sample = sg_sample_F.Children;
    sg_sample.YAxis.TickValues = [-2 -1 0]*1e-4;
    sg_sample.YAxis.TickLabels = [-20 -10 0];
    xlabel('x-x_{tip} [m]');
    ylabel('\Delta\epsilon_{xx} [10^{-3}]');
    xlim([-0.05 0.05]);
    pbaspect([8,5,1])
    legend('show', 'Location','northwest')
    
    %% A/A0
    Vx_and_A_F = figure;
    
    %-- superimpose Contact area
    x_minXtip_4spatial = PhediStructCell{EvNum}.BigPicRotStruct.x-PhediStructCell{EvNum}.PhediData.PhotoLocation;
    frontTime = PhediStructCell{EvNum}.BigPicRotStruct.frontTime_interp./PhediStructCell{EvNum}.BigPicRotStruct.fps;
    [~,I1] = min(abs(x_minXtip_4spatial));
    t_minTtip_4tmprl = PhediStructCell{EvNum}.BigPicRotStruct.t-frontTime(I1);
    [~,I2] = min(abs(t_minTtip_4tmprl));
    A_spatial = PhediStructCell{EvNum}.BigPicRotStruct.DataMatNorm(I2,:);
    A_temporal = PhediStructCell{EvNum}.BigPicRotStruct.DataMatNorm(:,I1);
    x_minXtip_4tmprl = signal_ChangeTime2Space_mapping(t_minTtip_4tmprl,PhediStructCell{EvNum}.PhediData.PhotoLocation,PhediStructCell{EvNum}.BigPicRotStruct);
    
    hold on;
    plot(x_minXtip_4tmprl,A_temporal,'-','Color',rgb('DeepSkyBlue'),'DisplayName','pix over time');
    % plot(x_minXtip_4spatial,A_spatial,'.-','Color',rgb('DodgerBlue'),'DisplayName','frame over space');
    ylabel('A/A_0');
    for l = 1:length(shift_vals)
        plot([1,1]*shift_vals(l),[0.78 1.05],'--','LineWidth',2, 'Color', shiftColor(l,:), 'DisplayName', sprintf('tip@shift %1.3f', num2str(shift_vals(l))));
    end
    hold off;
    
    xlim([-0.035 0.035]);
    ylim([0.78 1.05]*1e-0);
    pbaspect([8,5,1])
    Vx_and_A_F.Name = 'Vx_and_A';
    
end