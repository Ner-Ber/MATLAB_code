%% preferences
DoPlot = 1;

%% load
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
for iEv = 1:length(relevantEvents)
    PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
    PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
end

% relevantEvents=11;

%% prepare for saving
Cf_vec = nan(size(relevantEvents));

shifts_mean = nan(size(relevantEvents));
shifts_var = nan(size(relevantEvents));
shifts_skwns = nan(size(relevantEvents));
shifts_RMS_mean = nan(size(relevantEvents));
shifts_RMS_var = nan(size(relevantEvents));

strchs_mean = nan(size(relevantEvents));
strchs_var = nan(size(relevantEvents));
strchs_skwns = nan(size(relevantEvents));
strchs_RMS_mean = nan(size(relevantEvents));
strchs_RMS_var = nan(size(relevantEvents));

%% iterate on events
for N_i = 1:length(relevantEvents(:))
    
    n_prsntPhds = 5;
    
    if DoPlot
        figure;
        axes_shift = subplot(2,3,1);
        axes_stretch = subplot(2,3,2);
        axes_info = subplot(2,3,3);
        axes_ShftRes = subplot(2,3,4);
        axes_StrchRes = subplot(2,3,5);
        axes_RMShist = subplot(2,3,6);
        
        
        %% plot shifts
        [~, RMS_vec_shift, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStructCell{relevantEvents(N_i)},[],1,[],[],axes_shift);
        axis([[-0.07 0.01]*1 [-1 12]*1e-6]);
        %     pbaspect([5 2 1]);
        title('shifts');
        
        
        %% plot streches
        [~, RMS_vec_strech, stretch_vec] = phedi_UxBestFitVerticalStretch(PhediStructCell{relevantEvents(N_i)},[],1,[],[],axes_stretch);
        axis([[-0.07 0.01]*1 [-1 12]*1e-6]);
        %     pbaspect([5 2 1]);
        title('stretches');
        
        
        %% find which phedis not to present
        [B,I] = sort(shifts_vec);
        dn = (B(end)-B(1))/(n_prsntPhds-1);
        locs = B(1)+(1:(n_prsntPhds-2))*dn;
        diff_mat = bsxfun(@minus, B(:), locs(:)');
        [~,minIdxs] = min(abs(diff_mat));
        g = get(axes_shift,'children');
        f = get(axes_stretch,'children');
        presentIdx = [I(1);I(minIdxs);I(end)];
        presentIdxFull = [2*presentIdx(:)-1,2*presentIdx(:)];
        presentIdxFull = presentIdxFull(:);
        logicPlot = ~~(ones(length(g),1));
        for i=1:length(g)
            if sum(presentIdxFull==i)==0
                g(i).Visible='off';
                f(i).Visible='off';
            end
        end
        
        %% plot histogram
        axes(axes_RMShist); hold on;
        H_sh = histogram(RMS_vec_shift);
        H_st = histogram(RMS_vec_strech);
        M = max([H_sh.Values(:);H_st.Values(:)]);
        plot(mean(RMS_vec_shift).*[1 1],[0,M],'r');
        plot(mean(RMS_vec_strech).*[1 1],[0,M],'b');
        legend('shifts RMS','stretch RMS');
        
        
        if mean(RMS_vec_shift)>mean(RMS_vec_strech)
            str = ('stretch wins');
        elseif mean(RMS_vec_shift)<mean(RMS_vec_strech)
            str = ('shift wins');
        else
            str = ('TEKO!');
        end
        disp(str);
        
        
        %% write data
        axes(axes_info);
        text(0,0,PhediStructCell{relevantEvents(N_i)}.BigPicRotStruct.details,'FontSize',15);
        text(0,-0.5,['C_f=',num2str(PhediStructCell{relevantEvents(N_i)}.PhediData.Cf)],'FontSize',15);
        text(0,-1,[str,' by ',num2str(abs(mean(RMS_vec_shift)-mean(RMS_vec_strech)))],'FontSize',15);
        axis([-0.2 5 -1.5 0.5]);
        
        %% plot histograms of results
        axes(axes_ShftRes);
        histogram(shifts_vec);
        title(['\mu=',num2str(mean(shifts_vec)),'   \sigma^2=',num2str(var(shifts_vec)),'   skw=',num2str(skewness(shifts_vec))])
        axes(axes_StrchRes);
        histogram(stretch_vec,'FaceColor','b');
        title(['\mu=',num2str(mean(stretch_vec)),'   \sigma^2=',num2str(var(stretch_vec)),'   skw=',num2str(skewness(stretch_vec))])
        
    else
        [~, RMS_vec_shift, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStructCell{relevantEvents(N_i)},[],0);
        [~, RMS_vec_strech, stretch_vec] = phedi_UxBestFitVerticalStretch(PhediStructCell{relevantEvents(N_i)},[],0);
    end
    %% save data from loop
    Cf_vec(N_i) = PhediStructCell{relevantEvents(N_i)}.PhediData.Cf;
    
    shifts_mean(N_i) = mean(shifts_vec);
    shifts_var(N_i) = var(shifts_vec);
    shifts_skwns(N_i) = skewness(shifts_vec);
    shifts_RMS_mean(N_i) = mean(RMS_vec_strech);
    shifts_RMS_var(N_i) = var(RMS_vec_strech);
    
    strchs_mean(N_i) = mean(stretch_vec);
    strchs_var(N_i) = var(stretch_vec);
    strchs_skwns(N_i) = skewness(stretch_vec);
    strchs_RMS_mean(N_i) = mean(RMS_vec_shift);
    strchs_RMS_var(N_i) = var(RMS_vec_shift);
    %% save
    if DoPlot
        %         saveCurrentFigure('G:\Frics\2018-10-10\compareStretchAndShiftForDelta_i\compareAll_differential_FittingRegion',PhediStructCell{relevantEvents(N_i)}.BigPicRotStruct.details);
    end
end