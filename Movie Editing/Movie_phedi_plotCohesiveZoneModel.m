function ChsvAndAvgStrct = Movie_phedi_plotCohesiveZoneModel(DataStruct,varargin)
    
    %% set defaults
    %-- default cohesive model
    cohsvMdlDef = @(x) (1 - x.^2).*(0.5*tanh(3*x + 1) + 0.5) + (x + 1).^2.*(0.5*tanh(-3*x - 1) + 0.5);
    
    [cohsvMdl, plotSlopes, Handle] = setDefaults4function(varargin,...
        cohsvMdlDef, '11', []);
    
    plotDamaged = str2double(plotSlopes(end));
    slope2Plot = str2double(plotSlopes(1:end-1));
    
    %% general stuff
    SgData = DataStruct.SgData;
    PhediData = DataStruct.PhediData;
    
    scratchRegionMeters = DataStruct.Param.scratchRegionMeters;
    relevantSGs = find(scratchRegionMeters(1)<SgData.x_sg*1e-3 & SgData.x_sg*1e-3<scratchRegionMeters(2));
    if isfield(PhediData,'Cf')
        Cf = PhediData.Cf;
    else
        error('no Cf in data');
    end
    
    PhediSlopes = PhediData.slopeIncline;
    PhediSlopes = mean(PhediSlopes,1);
    PhediSlopes(~~mod(PhediSlopes,1)) = 0;
    Event = DataStruct.ExperimentData.eventNum;
    %--phedis to plot by number, pick relevant phedis index
    phdIndx = [];
    if plotDamaged
        phdIndx = cat(2,phdIndx,find(PhediSlopes==0));
    end
    if slope2Plot==-1
        phdIndx = cat(2,phdIndx,find(PhediSlopes==-1));
    elseif slope2Plot==1
        phdIndx = cat(2,phdIndx,find(PhediSlopes==1));
    else
        phdIndx = cat(2,phdIndx,find(PhediSlopes~=0));
    end
    phdIndx = unique(phdIndx);
    
    
    %% define handle of plot
    if isempty(Handle)
        figure('Name',['Event ',num2str(Event),' Avg. phedi and cohesive model']);
    else
        if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
            H = figure(Handle);
            H.Name = ['Event ',num2str(Event),' Avg. phedi and cohesive model'];
        elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
            axes(Handle);
        else
            error('''Handle'' isn''t axes nor figure');
        end
    end
    hold on;
    %% create models and average structure
    %-- create average phedi data
    AvgPhediStruct = phedi_averagePhedi(DataStruct,phdIndx);
    %-- find cohesive length
    %     [L_chosen, Dx] = phedi_fitCohesiveModel2Vel(PhediData.Cf,DataStruct.solAtSG.Gamma,AvgPhediStruct.Vel_mean_x,AvgPhediStruct.X_mean,cohsvMdl);
    Dx = 0; %0.005;
    L_chosen = 13e-3;
    %-- cohesive model on sg height
    CohesiveModelStructSg = CrackSolutionGeneralCohesive_InsertVariables_Neri(PhediData.Cf/DataStruct.Param.Cr,DataStruct.SgData.y_sg(8)*1e-3,...
        DataStruct.solAtSG.Gamma,'L',L_chosen,cohsvMdl);
    %--create model for interface
    CohesiveModelStructIntrfc = CrackSolutionGeneralCohesive_InsertVariables_Neri(PhediData.Cf/DataStruct.Param.Cr,1e-7,DataStruct.solAtSG.Gamma,...
        'L',L_chosen,cohsvMdl);
    
    %--- plot particle velocity from sg
    Uxx_i = (SgData.Uxx(:,8)-mean(SgData.Uxx(1:100,8)))*1e-3;
    SG_x_mins_x_tip = SgData.x_mins_x_tips(:,8);
    plot(SG_x_mins_x_tip,-Cf.*Uxx_i,...
        '.-','Color',rgb('DarkGreen'),'MarkerSize',10,'LineWidth',1,'DisplayName','sg particle vel.');
    %-- plot particle velocity LEFM
    plot(DataStruct.solAtSG.x,-Cf.*DataStruct.solAtSG.Uxx,...
        '.-','Color',rgb('LightSeaGreen'),'LineWidth',2,'DisplayName','LEFM @sg');
    %-- plot cohsv on interface
    plot(CohesiveModelStructSg.x,CohesiveModelStructSg.vx,'Color',rgb('steelBlue'),'LineWidth',3,'DisplayName','Cohesive Model @sg');
    
    
    %--plot shaded variance
    XX = [AvgPhediStruct.X_mean;flipud(AvgPhediStruct.X_mean)];
    YY = [AvgPhediStruct.Vel_mean_x-sqrt(AvgPhediStruct.Vel_var_x);flipud(AvgPhediStruct.Vel_mean_x+sqrt(AvgPhediStruct.Vel_var_x))];
    %-- eliminate Nan
    YY = YY(~isnan(XX));
    XX = XX(~isnan(XX));
    XX = XX(~isnan(YY));
    YY = YY(~isnan(YY));
    fill(XX-Dx,YY,rgb('PeachPuff'),'linestyle','none','FaceAlpha',0.5,'DisplayName','phedi Vel std');
    %-- plot mean phedi
    plot(AvgPhediStruct.X_mean-Dx,AvgPhediStruct.Vel_mean_x,'Color',rgb('Crimson'),'LineWidth',3,'DisplayName','mean phedi');
    %-- plot cohsv on interface
    plot(CohesiveModelStructIntrfc.x,CohesiveModelStructIntrfc.vx,'Color',rgb('mediumOrchid'),'LineWidth',3,'DisplayName','Cohesive Model @intrface');
    %-- plot LEFM on interface
    plot(DataStruct.solAtInter.x,DataStruct.solAtInter.vx,'Color',rgb('silver'),'LineWidth',2,'DisplayName','LEFM @intrface');
    
    
    
    
    
    
    %--- make plot pretty
    xlim([-0.05 0.05]);
    ylim([-0.02 1.3*max(YY(abs(XX)<=0.05))])
    xlabel('x-x_{tip} [m]');
    ylabel('velocity [m/s]');
    title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)])
    hold off;
    
    %% save to structure
    if nargout>0
        ChsvAndAvgStrct.AvgPhediStruct = AvgPhediStruct;
        ChsvAndAvgStrct.CohesiveModelStructSg = CohesiveModelStructSg;
        ChsvAndAvgStrct.CohesiveModelStructIntrfc = CohesiveModelStructIntrfc;
    end
end