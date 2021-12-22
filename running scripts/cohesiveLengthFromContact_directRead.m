[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
D = dir_names;
Experiments = find(~cellfun(@isempty,strfind(D,'-')));
% exper = my_dir;
% Experiments = 1:length(exper);
% findLocs =[0.02 0.03 0.04 0.06 0.07 0.08 0.1 0.11 0.13 0.14];
% findLocs =[0.04 0.06 0.07 0.08 0.1 0.11 0.13];
findLocs = linspace(0.02,0.14,30);
Cols4Locs = distinguishable_colors(length(findLocs));
Marker4Exp = {'*','o','s','d','^','v','x','p','h','+','>','<','*','o','s','d','^','v','x','p','h','+','>','<','*','o','s','d','^','v','x','p','h','+','>','<'};
Size4Exp = [28*ones(1,12) 36*ones(1,12) 44*ones(1,12)];

XcFig = figure; hold on;
plot([1 1]*Cr,[0 4e-3],'--','Color',[0.2 0.2 0.2],'LineWidth',2);
text(Cr+10,1.5e-3,'C_R','Color',[0.2 0.2 0.2]);
plot([1 1]*Cs,[0 4e-3],'--','Color',[0.8 0.8 0.8],'LineWidth',2);
text(Cs+10,1.5e-3,'C_S','Color',[0.8 0.8 0.8]);
colormap(Cols4Locs)
Ch = colorbar;
Ch.Ticks = movmean(Ch.Ticks,2,'Endpoints','discard');
TL = cellfun(@num2str,num2cell(findLocs),'UniformOutput',0);
Ch.TickLabels = TL;
Ch.Label.String = 'location on block [m]';

counter = 1;
clear allChsvS;
for iExp = 1:length(Experiments)
    PhData = phantomReadMeta(fullfile(D{iExp},'Ph'));
%     PhData = phantomReadMeta(fullfile(D{iExp},'PhBig'));
    relevantEvents = 1:PhData.NumEvents;
    for I=1:length(relevantEvents)
        BigPicRotStruct = phantom_SmallPicStruct(exper{iExp}, I,[],[],[],[],[],4:5);
%         BigPicRotStruct = phantom_BigPicStruct(D{iExp}, I,[],[],[],[],[],4:5);
        if phantom_isPrecursor(BigPicRotStruct);
            continue
        end
        Cf_vec = nan(length(findLocs),1);
        Xc_vec = nan(length(findLocs),1);
%         Xc_vec_alter = nan(length(findLocs),1);
        for f = 1:length(findLocs)
            try
            [~,iLoc] = min(abs(smooth(BigPicRotStruct.frontVelLoc_interpM,5)-findLocs(f)));
            Cf_vec(f) = BigPicRotStruct.frontVel_interpMperS(iLoc);
%             ChsvFromCntct = ROT_findXcFromContact_fromPix(BigPicRotStruct,findLocs(f));
            ChsvFromCntct = ROT_findXcFromContact_fromFrame(BigPicRotStruct,findLocs(f));
            Xc_vec(f) = ChsvFromCntct.Xc;
%             Xc_vec_alter(f) = -ChsvFromCntct.Tc*Cf_vec(f);
            catch
                warning(['loc ',num2str(findLocs(f)),', Ev ',num2str(relevantEvents(I)),' in exp ',D{Experiments(iExp)},' didn''t work']);
            end
            allChsvS(counter) = ChsvFromCntct;
            counter = counter+1;
        end
        Name = [BigPicRotStruct.details];
        figure(XcFig);
        scatter(Cf_vec,abs(Xc_vec),Size4Exp(iExp),Cols4Locs,Marker4Exp{iExp},'DisplayName',Name,'LineWidth',2);
%         plot(Cf_vec,abs(Xc_vec_alter),'.','Color',[0.2 0.2 0.2]);
    end
    
    
end

