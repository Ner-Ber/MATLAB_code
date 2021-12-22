function Movie_phedi_plotWithSg(DataStruct,varargin)
    % Movie_phedi_plotWithSg(DataStruct,plotSlopes, plotNames, reduceSin, originlTime, plotInBlock, axesHandle)
    %
    %Movie_phedi_plotWithSg will plot phedi tracking given the output of the
    %function 'Movie_phedi_from_folder_2_data';
    %
    % VARARGIN = DAFAULTS:
    % plotSlopes='11'       (two digits:
    %                       1st-indicates slope to plot (+/-1) 0 means all slopes.
    %                       2nd-logical indicating if plot damaged.
    % 
    % plotNames={'ROT','loc','vel','phase','locSpace','velSpace','phediUxx','avgPhediLoc','avgPhediVel'}
    %                       This is a cell array containing names of plot you want.
    %                       relevant only in 'sep' mode.
    % 
    % reduceSin=1           reduce a fitted sin form from the phedi
    %                       displacements
    %
    % originlTime=0         false - use shifted to t-t_tip vector.
    %                       true - use original time vector.
    %
    % plotInBlock=1         positive - uses logicals to plot only parts of
    %                       vectors that are NOT extrapolated.
    %
    % Handle=[]             notempty - plot in the given axes. (works only
    %                       in a single plotting option)
    
    %% set defaults
    [plotSlopes, plotNames, reduceSin, originlTime, plotInBlock, Handle] = setDefaults4function(varargin,...
        '11', {'ROT','loc','vel','phase','locSpace','velSpace','phediUxx'},1, 0, 1, []);
    
    
    plotDamaged = str2double(plotSlopes(end));
    slope2Plot = str2double(plotSlopes(1:end-1));
    
    
    %% general stuff
    [Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
    AsperityData = DataStruct.AsperityData;
    SgData = DataStruct.SgData;
    PhediData = DataStruct.PhediData;
    BigPicRotStruct = DataStruct.BigPicRotStruct; 
    
    scratchRegionMeters = DataStruct.Param.scratchRegionMeters;
    relevantSGs = find(scratchRegionMeters(1)<SgData.x_sg*1e-3 & SgData.x_sg*1e-3<scratchRegionMeters(2));
    SgColors = MyVaryColor(size(SgData.Uxx,2));
    if isfield(PhediData,'Cf')
        Cf = PhediData.Cf;
    else
        error('no Cf in data');
    end
    
    PhediLocation = PhediData.PhediLocation;
    PhediSlopes = PhediData.slopeIncline;
            PhediSlopes = mean(PhediSlopes,1);
    PhediSlopes(~~mod(PhediSlopes,1)) = 0;
    Event = DataStruct.ExperimentData.eventNum;
    center_of_mass = DataStruct.AsperityData.center_of_mass;
    res = DataStruct.ExperimentData.res;
    CM_time = DataStruct.AsperityData.timeCount;
    Nphd = size(PhediLocation,2);       % total number of phedis
    %---- manage plot properties:
    %--colors:
    FigColors = MyVaryColor(Nphd);
    FigColors(PhediSlopes==0,:) = repmat([1 0 0],nnz(PhediSlopes==0),1);    % set damages to red
    %--names:
    LGND = cellfun(@num2str,num2cell(PhediData.measuredPhedisFromPlot),'UniformOutput',0);
    %--markers:
    allMarks = {'*-','.-','--'};
    marksLookup = PhediSlopes;
    marksLookup(PhediSlopes==-1) = 2;
    marksLookup(PhediSlopes==0) = 3;
    FigMarks = allMarks(marksLookup);
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
    
    
    
    
    %% plot the heat map with phedis
    if any(strcmp(plotNames,'ROT'))
        if isempty(Handle)
            figure('Name',['Event ',num2str(Event),' RowOver w/ Phedis']);
        else
            if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
                H = figure(Handle);
                H.Name = ['Event ',num2str(Event),' RowOver w/ Phedis'];
            elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
                axes(Handle);
            else
                error('''Handle'' isn''t axes nor figure');
            end
        end
        hold on;
        shiftByLocation = PhediData.PhotoLocation-mean(AsperityData.spatialVec);
        xImagesc = [AsperityData.spatialVec(1) AsperityData.spatialVec(end)]+shiftByLocation;
        yImagesc = [AsperityData.timeCount(1), AsperityData.timeCount(end)];
        imagesc(xImagesc,yImagesc,AsperityData.RowOverTime);
        
        %--- plot it:
        for i = phdIndx(:)'
%                 plot(PhediLocation(:,i)+PhediData.measuredPhedisFromPlot(i)+shiftByLocation,...
%                     PhediData.timeVec,...
%                     FigMarks{i},'Color',FigColors(i,:),'DisplayName',LGND{i});
                
            plot(PhediLocation(:,i)+shiftByLocation,...
                PhediData.timeVec,...
                FigMarks{i},'Color',FigColors(i,:),'DisplayName',LGND{i});
        end
        
        
        ylabel('tims [s]');
        xlabel('approx. location on block [m]');
        title({'row over time',[DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]});
        xlim(xImagesc); ylim(yImagesc);

        %--- add center of mass trajectories
        plot(center_of_mass/res+shiftByLocation,CM_time);
        
    end
    
    %% plot phedi location
    if any(strcmp(plotNames,'loc'))
        Phedis2Zero = bsxfun(@minus,PhediLocation,mean(PhediLocation(1:30,:),'omitnan'));
        if reduceSin
            PhediLocReduced = phedi_reduceSinFromLocation(PhediData);
            PhediLocReduced = bsxfun(@minus,PhediLocReduced,mean(PhediLocReduced(1:30,:),'omitnan'));
        else
            PhediLocReduced = Phedis2Zero;
        end
        
        if isempty(Handle)
            figure('Name',['Event ',num2str(Event),' RowOver w/ Phedis']);
        else
            if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
                H = figure(Handle);
                H.Name = ['Event ',num2str(Event),' Phedis loc w/ Uxx'];
            elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
                axes(Handle);
            else
                error('''Handle'' isn''t axes nor figure');
            end
        end
        hold on;
        
        if originlTime
            timeVec = PhediData.timeVec;
            timeVec = repmat(timeVec(:),1,Nphd);
        else
            timeVec = PhediData.t_mins_t_tip;
        end
        
        %--- plot it:
        for i = phdIndx(:)'
                plot(timeVec(:,i),PhediLocReduced(:,i),...
                    FigMarks{i},'Color',FigColors(i,:),'DisplayName',LGND{i});
        end
        
        %--- add cm location
        %     center_of_mass_zero = bsxfun(@minus,center_of_mass,mean(center_of_mass(10:100,:),1));
        %     hold on;
        %     plot(CM_time,center_of_mass_zero/res,'LineWidth',1.5);
        
        xlabel('tims [s]');
        ylabel('location of phedi [m]');
        title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]);
        xlim([timeVec(1), timeVec(end)]);
        hold off;
    end
    
    
    %% plot phedi location as function of location
    if any(strcmp(plotNames,'locSpace'))
        Phedis2Zero = bsxfun(@minus,PhediLocation,mean(PhediLocation(1:30,:),'omitnan'));
        if reduceSin
            PhediLocReduced = phedi_reduceSinFromLocation(PhediData);
            PhediLocReduced = bsxfun(@minus,PhediLocReduced,mean(PhediLocReduced(1:30,:),'omitnan'));
        else
            PhediLocReduced = Phedis2Zero;
        end
        
        if isempty(Handle)
            figure('Name',['Event ',num2str(Event),' RowOver w/ Phedis']);
        else
            if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
                H = figure(Handle);
                H.Name = ['Event ',num2str(Event),' Phedis loc Vs. location'];
            elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
                axes(Handle);
            else
                error('''Handle'' isn''t axes nor figure');
            end
        end
        hold on;
        
        %--- set spatial axis
        spaceAxis = PhediData.x_mins_x_tip;
        
        %-- set limits for relevant plotting
        if plotInBlock && isfield(PhediData,'inBlockLogicals')
            inBlockLogicals = logical(PhediData.inBlockLogicals);
        else
            inBlockLogicals = true(size(spaceAxis));
        end
        
        %--- plot it:
        Errs = repmat(1.0811e-06,size(inBlockLogicals,1),1);
        for i = phdIndx(:)'
            plot(spaceAxis(inBlockLogicals(:,i),i),PhediLocReduced(inBlockLogicals(:,i),i),...
                FigMarks{i},'Color',FigColors(i,:),'DisplayName',LGND{i});
            
%             errorbar(spaceAxis(inBlockLogicals(:,i),i),PhediLocReduced(inBlockLogicals(:,i),i),Errs(inBlockLogicals(:,i))/2,...
%                 FigMarks{i},'Color',FigColors(i,:),'DisplayName',LGND{i});
        end
        
        xlabel('x-x_{tip} [m]');
        ylabel('location of phedi [m]');
        title('location of phedis');
        
        title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]);
        xlim([-0.05 0.05]);
        hold off;
    end
    %% plot speed from phedis and sg
    if any(strcmp(plotNames,'vel'))
        if isempty(Handle)
            figure('Name',['Event ',num2str(Event),' RowOver w/ Phedis']);
        else
            if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
                H = figure(Handle);
                H.Name = ['Event ',num2str(Event),' velocity'];
            elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
                axes(Handle);
            else
                error('''Handle'' isn''t axes nor figure');
            end
        end
        hold on;
        
        %--- set time axis
        if originlTime
            [PhediVelocity, timeVec_4vel] = deal(PhediData.PhediVelocity,PhediData.timeVec_4vel);
            timeVec_4vel = repmat(timeVec_4vel(:),1,size(PhediVelocity,2));
        else
            [PhediVelocity, timeVec_4vel] = deal(PhediData.PhediVelocity,PhediData.t_mins_t_tip_4vel);
        end
        
        %--- plot it:
        for i = phdIndx(:)'
            plot(timeVec_4vel(:,i),PhediVelocity(:,i),...
                FigMarks{i},'Color',FigColors(i,:),'DisplayName',LGND{i});
        end        

        xlabel('t-t_{tip} [s]');
        ylabel('velocity of phedi [m/s]');
        ylim([-2, 2]);
        
        %--- plot particle velocity from sg
        for i = relevantSGs
            %--- particle velocity from sg:
            Uxx_i = (SgData.Uxx(:,8)-mean(SgData.Uxx(1:100,8)))*1e-3;
            CurrentParticleVelocity = sg_calcParticleVelDiffCf(SgData.t_mins_t_tips(:,i), SgData.x_sg(i)*1e-3, DataStruct.BigPicRotStruct,Uxx_i);
            %--- find t_tip:
            if originlTime
                SG_timeVec = (SgData.t./1000);
            else
                SG_timeVec = SgData.t_mins_t_tips(:,i);
            end
%             plot(SG_timeVec,CurrentParticleVelocity,...
%                 '.-','Color',rgb('DarkGreen'),'MarkerSize',10,'LineWidth',1);
            plot(SG_timeVec,-Cf.*Uxx_i,...
                '.-','Color',rgb('Orange'),'MarkerSize',10,'LineWidth',1);
        end
        ylabel('particle velocity (m/s)');
        legend(LGND,'Location','eastoutside');
        legend('off');
        xlim([-1 1]*2e-4);
        hold off;
        title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]);
    end
    %% plot speed as function of location (phase space plot)
    if any(strcmp(plotNames,'phase'))
        if isempty(Handle)
            figure('Name',['Event ',num2str(Event),' RowOver w/ Phedis']);
        else
            if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
                H = figure(Handle);
                H.Name = ['Event ',num2str(Event),' phase space'];
            elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
                axes(Handle);
            else
                error('''Handle'' isn''t axes nor figure');
            end
        end
        hold on;
        
        PhediVelocity = PhediData.PhediVelocity;
        Phedis2Zero = bsxfun(@minus,PhediLocation,mean(PhediLocation(1:30,:),'omitnan'));
        if reduceSin
            PhediLocReduced = phedi_reduceSinFromLocation(PhediData);
            PhediLocReduced = bsxfun(@minus,PhediLocReduced,mean(PhediLocReduced(1:30,:),'omitnan'));
        else
            PhediLocReduced = Phedis2Zero;
        end
        PhediLoc4Phase = movmean(PhediLocReduced,2,1,'Endpoints','discard');
        %--- plot it:
        for i = phdIndx(:)'
            plot(PhediLoc4Phase(:,i),PhediVelocity(:,i),...
                FigMarks{i},'Color',FigColors(i,:),'DisplayName',LGND{i});
        end 

        xlabel('relative location [m]');
        ylabel('velocity of phedi [m/s]');
        %     title('location=f(velocity)');
        title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)]);
        
    end
    
    %% plot particle velocity is space with cohesive zone model
    if any(strcmp(plotNames,'velSpace'))
        if isempty(Handle)
            figure('Name',['Event ',num2str(Event),' RowOver w/ Phedis']);
        else
            if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
                H = figure(Handle);
                H.Name = ['Event ',num2str(Event),' vel in space+ cohesive'];
            elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
                axes(Handle);
            else
                error('''Handle'' isn''t axes nor figure');
            end
        end
        hold on;
        [PhediVelocity, spaceVec_4vel] = deal(PhediData.PhediVelocity,PhediData.x_mins_x_tip_4vel);
        
        %-- set limits for relevant plotting
        if plotInBlock && isfield(PhediData,'inBlockLogicals') && size(PhediVelocity,2)==size(PhediData.x_mins_x_tip,2)
            inBlockLogicals = logical(PhediData.inBlockLogicals);
        else
            inBlockLogicals = true(size(PhediVelocity));
        end
        
        %--- plot it:
        for i = phdIndx(:)'
            plot(spaceVec_4vel(inBlockLogicals(:,i),i),PhediVelocity(inBlockLogicals(:,i),i),...
                FigMarks{i},'Color',FigColors(i,:),'DisplayName',LGND{i});
        end 
        
        
        %--- plot particle velocity from sg
        for i = relevantSGs
            %--- particle velocity from sg:
            Uxx_i = (SgData.Uxx(:,8)-mean(SgData.Uxx(1:100,8)))*1e-3;
            CurrentParticleVelocity = sg_calcParticleVelDiffCf(SgData.t_mins_t_tips(:,i), SgData.x_sg(i)*1e-3, DataStruct.BigPicRotStruct,Uxx_i);
            SG_x_mins_x_tip = SgData.x_mins_x_tips(:,i);
%             plot(SG_x_mins_x_tip,CurrentParticleVelocity,...
%                 '.-','Color',rgb('DarkGreen'),'MarkerSize',10,'LineWidth',1,'DisplayName','sg particle vel.');
            plot(SG_x_mins_x_tip,-Cf.*Uxx_i,...
                '.-','Color',rgb('DarkGreen'),'MarkerSize',10,'LineWidth',1,'DisplayName','sg particle vel.');
        end
        
        %--- plot particle velocity at interface from cohesive zone model
%         if isfield(DataStruct,'CohesiveModelStruct')
%             sol=DataStruct.CohesiveModelStruct;
%             x_4_Ux_dot = movmean(sol.x,2,'Endpoints','discard');
%             t = -sol.x/Cf;
%             Ux_dot = diff(sol.Ux)./diff(t);
%             plot(x_4_Ux_dot,Ux_dot,'Color',[1 0 1],'LineWidth',2,'DisplayName','cohesive model at interface');
%         else
%             plot(nan,nan,'Color',[1 0 1],'LineWidth',2);
%         end
        %--- make plot pretty
        xlabel('x-x_{tip} [m]');
        ylabel('velocity [m/s]');
        title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)])
        hold off;
        
        
    end
    
    %% plot strain calculated from phedis
    if any(strcmp(plotNames,'phediUxx'))
        phedisLoc4Uxx = bsxfun(@plus, PhediData.PhediLocation,PhediData.measuredPhedisFromPlot');
        % phedisNums = [4 6 8 10 12 14 16 18 22 24 26 30];
        %     phedisSlopesFiction = zeros(1,size(PhediData.PhediLocation,2));
        % phedisSlopesFiction(phedisNums) = 1;
        %     phedisSlopesFiction(1:length(phedisSlopesFiction)) = 1;
        
        %     [UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx(phedisLoc4Uxx,[],[],0.15e-3,phedisSlopesFiction,1);
        [UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx_ofPairs(phedisLoc4Uxx);
        
        locations4Legend = cellfun(@num2str,num2cell(round(initial_locations,5)),'UniformOutput',0);
        
        figure;
        pairsLGND = {};
        hold on;
        for ii = 1:size(phediPairs,2)
            plot(PhediData.timeVec,UxxFromPhedi(:,ii),'Linewidth',1.5);
            pairsLGND{ii} = [locations4Legend{1,ii},' & ' locations4Legend{2,ii}];
        end
        UxxFromPhediMEAN = mean(UxxFromPhedi,2);
        plot(PhediData.timeVec,UxxFromPhediMEAN,'k','Linewidth',1.5);
        pairsLGND{ii+1} = 'mean Uxx';
        legend(pairsLGND,'Location','WestOutside');
        title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf)])
        hold off;
    end
   
   %% plot average phedi location
    if any(strcmpi(plotNames,'avgPhediLoc'))
        if isempty(Handle)
            figure('Name',['Event ',num2str(Event),' Avg. phedi location']);
        else
            if isa(Handle,'matlab.ui.Figure')   % if handle is a figure
                H = figure(Handle);
                H.Name = ['Event ',num2str(Event),' Avg. phedi location'];
            elseif isa(Handle,'matlab.graphics.axis.Axes')   % if handle is a axes
                axes(Handle);
            else
                error('''Handle'' isn''t axes nor figure');
            end
        end
        
        %-- create average phedi data
        AvgPhediStruct = phedi_averagePhedi(DataStruct,phdIndx);
%         yyaxis left;
        hold on;
        
        
        %--- plot particle velocity from sg
        Uxx_i = (SgData.Uxx(:,8)-mean(SgData.Uxx(1:100,8)))*1e-3;
        SG_x_mins_x_tip = SgData.x_mins_x_tips(:,8);
        
        Ux_tmp=-cumsum(Uxx_i.*gradient(SG_x_mins_x_tip));% [m] for slip multiply by 2.
        [~,ii] = min(abs(SG_x_mins_x_tip));
        Ux_sg = -(Ux_tmp-Ux_tmp(ii));
%         Ux_sg(SG_x_mins_x_tip>=0) = 0;
        solAtSG = CrackSolutionForh_GammaChange(DataStruct.PhediData.Cf/Cr,-0.5,3.5e-3,DataStruct.solAtSG.Gamma,(-0.5:0.001:0.5));
        %-- zero both displacement functions
        zeroPoint = min(max(SG_x_mins_x_tip),max(solAtSG.x));
        [~,I] = min(abs(SG_x_mins_x_tip-zeroPoint));
        Ux_sg = Ux_sg-Ux_sg(I);
        
        [~,I] = min(abs(solAtSG.x-zeroPoint));
        UxSol = solAtSG.Ux-solAtSG.Ux(I);
        
        plot(SG_x_mins_x_tip,Ux_sg,...
            '.-','Color',rgb('DarkGreen'),'MarkerSize',10,'LineWidth',1,'DisplayName','sg particle vel.');
        %-- plot particle velocity LEFM
        plot(solAtSG.x,UxSol,...
            '.-','Color',rgb('LightSeaGreen'),'LineWidth',2,'DisplayName','LEFM @sg');

        
        %--plot shaded variance
        XX = [AvgPhediStruct.X_mean;flipud(AvgPhediStruct.X_mean)];
%         YY = [AvgPhediStruct.Vel_mean_x-sqrt(AvgPhediStruct.Vel_var_x);flipud(AvgPhediStruct.Vel_mean_x+sqrt(AvgPhediStruct.Vel_var_x))];
        YY = [AvgPhediStruct.Loc_mean_x-sqrt(AvgPhediStruct.Loc_var_x);flipud(AvgPhediStruct.Loc_mean_x+sqrt(AvgPhediStruct.Loc_var_x))];
        %-- eliminate Nan
        YY = YY(~isnan(XX));
        XX = XX(~isnan(XX));
        XX = XX(~isnan(YY));
        YY = YY(~isnan(YY));
        fill(XX,YY,rgb('PeachPuff'),'linestyle','none','FaceAlpha',0.5,'DisplayName','phedi Vel std');
        %-- plot mean phedi
        plot(AvgPhediStruct.X_mean,AvgPhediStruct.Loc_mean_x,'-','Color',rgb('Crimson'),'LineWidth',3,'DisplayName','mean phedi');
        
        %-- plot LEFM on interface
        plot(DataStruct.solAtInter.x,DataStruct.solAtInter.Ux,'-','Color',rgb('silver'),'LineWidth',2,'DisplayName','LEFM @intrface');

        %--- make plot pretty
        xlim([-0.05 0.05]);
        ylim([-0.02*1e-6 1.3*max(YY(abs(XX)<=0.05))])
        xlabel('x-x_{tip} [m]');
        ylabel('u_x [m/s]');
        title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf),' loc'])
        hold off;
        
        
%         yyaxis right;
%         %-- superimpose Contact area
%         x_minXtip_4spatial = BigPicRotStruct.x-PhediData.PhotoLocation;
%         frontTime = BigPicRotStruct.frontTime_interp./BigPicRotStruct.fps;
%         [~,I1] = min(abs(x_minXtip_4spatial));
%         t_minTtip_4tmprl = BigPicRotStruct.t-frontTime(I1);
%         [~,I2] = min(abs(t_minTtip_4tmprl));
%         A_spatial = BigPicRotStruct.DataMatNorm(I2,:);
%         A_temporal = BigPicRotStruct.DataMatNorm(:,I1);
%         x_minXtip_4tmprl = signal_ChangeTime2Space_mapping(t_minTtip_4tmprl,PhediData.PhotoLocation,BigPicRotStruct);
%         
%         hold on;
%         plot(x_minXtip_4tmprl,A_temporal,'.-','Color',rgb('DeepSkyBlue'),'DisplayName','pix over time');
%         plot(x_minXtip_4spatial,A_spatial,'.-','Color',rgb('DodgerBlue'),'DisplayName','frame over space');
%         ylabel('A/A_0');
%         hold off;
        
        
        
    end
    
    %% plot average phedi velocity
    if any(strcmpi(plotNames,'avgPhediVel'))
        if isempty(Handle)
            figure('Name',['Event ',num2str(Event),' Avg. phedi velocity']);
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
        
        %-- create average phedi data
        AvgPhediStruct = phedi_averagePhedi(DataStruct,phdIndx);
%         yyaxis left;
        hold on;
        
        
        %--- plot particle velocity from sg
        Uxx_i = (SgData.Uxx(:,8)-mean(SgData.Uxx(1:100,8)))*1e-3;
        SG_x_mins_x_tip = SgData.x_mins_x_tips(:,8);
        plot(SG_x_mins_x_tip,-Cf.*Uxx_i,...
            '.-','Color',rgb('DarkGreen'),'MarkerSize',10,'LineWidth',1,'DisplayName','sg particle vel.');
        %-- plot particle velocity LEFM
        plot(DataStruct.solAtSG.x,-Cf.*DataStruct.solAtSG.Uxx,...
            '.-','Color',rgb('LightSeaGreen'),'LineWidth',2,'DisplayName','LEFM @sg');
%         plot(DataStruct.solAtSG.x,DataStruct.solAtSG.vx,...
%             '.-','Color',rgb('LightSeaGreen'),'LineWidth',2,'DisplayName','LEFM @sg');

        
        %--plot shaded variance
        XX = [AvgPhediStruct.X_mean;flipud(AvgPhediStruct.X_mean)];
        YY = [AvgPhediStruct.Vel_mean_x-sqrt(AvgPhediStruct.Vel_var_x);flipud(AvgPhediStruct.Vel_mean_x+sqrt(AvgPhediStruct.Vel_var_x))];
        %-- eliminate Nan
        YY = YY(~isnan(XX));
        XX = XX(~isnan(XX));
        XX = XX(~isnan(YY));
        YY = YY(~isnan(YY));
        fill(XX,YY,rgb('PeachPuff'),'linestyle','none','FaceAlpha',0.5,'DisplayName','phedi Vel std');
        %-- plot mean phedi
        plot(AvgPhediStruct.X_mean,AvgPhediStruct.Vel_mean_x,'-','Color',rgb('Crimson'),'LineWidth',3,'DisplayName','mean phedi');
        
        %-- plot LEFM on interface
        plot(DataStruct.solAtInter.x,DataStruct.solAtInter.vx,'-','Color',rgb('silver'),'LineWidth',2,'DisplayName','LEFM @intrface');

        %--- make plot pretty
        xlim([-0.05 0.05]);
        ylim([-0.02 1.3*max(YY(abs(XX)<=0.05))])
        xlabel('x-x_{tip} [m]');
        ylabel('velocity [m/s]');
        title([DataStruct.BigPicRotStruct.details, ' Cf=',num2str(PhediData.Cf),' vel'])
        hold off;
        
        
        yyaxis right;
        %-- superimpose Contact area
        x_minXtip_4spatial = BigPicRotStruct.x-PhediData.PhotoLocation;
        frontTime = BigPicRotStruct.frontTime_interp./BigPicRotStruct.fps;
        [~,I1] = min(abs(x_minXtip_4spatial));
        t_minTtip_4tmprl = BigPicRotStruct.t-frontTime(I1);
        [~,I2] = min(abs(t_minTtip_4tmprl));
        A_spatial = BigPicRotStruct.DataMatNorm(I2,:);
        A_temporal = BigPicRotStruct.DataMatNorm(:,I1);
        x_minXtip_4tmprl = signal_ChangeTime2Space_mapping(t_minTtip_4tmprl,PhediData.PhotoLocation,BigPicRotStruct);
        
        hold on;
        plot(x_minXtip_4tmprl,A_temporal,'.-','Color',rgb('DeepSkyBlue'),'DisplayName','pix over time');
        plot(x_minXtip_4spatial,A_spatial,'.-','Color',rgb('DodgerBlue'),'DisplayName','frame over space');
        ylabel('A/A_0');
        hold off;
        
        
    end
    
end