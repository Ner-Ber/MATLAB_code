[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
D = dir_names;
Experiments = find(~cellfun(@isempty,strfind(D,'.mat')));
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

% Experiments = Experiments(7:end);
Experiments = Experiments(end);
for iExp = 1:length(Experiments)
%     load(D{Experiments(iExp)});
    relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
    relevantEvents = relevantEvents(relevantEvents~=1);
%     relevantEvents=3;
    for I=1:length(relevantEvents)
        Cf_vec = nan(length(findLocs),1);
        Xc_vec = nan(length(findLocs),1);
        Xc_vec_alter = nan(length(findLocs),1);
        for f = 1:length(findLocs)
            try
            [~,iLoc] = min(abs(smooth(PhediStructCell{relevantEvents(I)}.BigPicRotStruct.frontVelLoc_interpM,5)-findLocs(f)));
            Cf_vec(f) = PhediStructCell{relevantEvents(I)}.BigPicRotStruct.frontVel_interpMperS(iLoc);
%             ChsvFromCntct = ROT_findXcFromContact_fromFrame(PhediStructCell{relevantEvents(I)}.BigPicRotStruct,findLocs(f));
            ChsvFromCntct = ROT_findXcFromContact_fromPix(PhediStructCell{relevantEvents(I)}.BigPicRotStruct,findLocs(f));
            Xc_vec(f) = ChsvFromCntct.Xc;
            Xc_vec_alter(f) = -ChsvFromCntct.Tc*Cf_vec(f);
            catch
                warning(['loc ',num2str(findLocs(f)),', Ev ',num2str(relevantEvents(I)),' in exp ',D{Experiments(iExp)},' didn''t work']);
            end
        end
        Name = [PhediStructCell{relevantEvents(I)}.BigPicRotStruct.details];
        figure(XcFig);
        scatter(Cf_vec,abs(Xc_vec),Size4Exp(iExp),Cols4Locs,Marker4Exp{iExp},'DisplayName',Name,'LineWidth',2);
        plot(Cf_vec,abs(Xc_vec_alter),'.','Color',[0.2 0.2 0.2]);
    end
    
    
end

