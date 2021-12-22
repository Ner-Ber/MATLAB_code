[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
D = dir_names;
Experiments = find(~cellfun(@isempty,strfind(D,'-')));
Experiments=1;
% exper = my_dir;
% Experiments = 1:length(exper);
% findLocs =[0.02 0.03 0.04 0.06 0.07 0.08 0.1 0.11 0.13 0.14];
% findLocs =[0.04 0.06 0.07 0.08 0.1 0.11 0.13];
% findLocs = linspace(0.02,0.14,30);
% Cols4Locs = distinguishable_colors(length(findLocs));
allChsvS = struct([]);
for iExp = 1:length(Experiments)
        PhData = phantomReadMeta(fullfile(D{iExp},'Ph'));
%     PhData = phantomReadMeta(fullfile(D{iExp},'PhBig'));
    relevantEvents = 2:PhData.NumEvents;
    for I=1:length(relevantEvents)
        try
            BigPicRotStruct = phantom_SmallPicStruct(exper{iExp}, I,[],[],[],[],[],4:5);
%             BigPicRotStruct = phantom_BigPicStruct(D{iExp}, relevantEvents(I),[],[],[],[],[],4:5);
            if phantom_isPrecursor(BigPicRotStruct);
                continue
            end
            ChsvFromCntct = ROT_findXcFromContact(BigPicRotStruct);
            ChsvFromCntct.CfVec = BigPicRotStruct.frontVel_interpMperS;
            ChsvFromCntct.Exp = D{iExp};
            ChsvFromCntct.Ev = relevantEvents(I);
            allChsvS = cat(1,allChsvS,ChsvFromCntct);
        catch
            disp('');
        end
    end
end


%% plot
Marker4Exp = {'*','o','s','d','^','v','x','p','h','+','>','<','*','o','s','d','^','v','x','p','h','+','>','<','*','o','s','d','^','v','x','p','h','+','>','<'};

all_Xc = arrayfun(@(S) movmean(S.Xc,2,'Endpoints','discard'), allChsvS,'UniformOutput',0);
all_Xc = cell2mat(all_Xc);
all_Cf = arrayfun(@(S) S.CfVec, allChsvS,'UniformOutput',0);
all_Cf = all_Cf';
all_Cf = cell2mat(all_Cf);
all_Xc =all_Xc';
Xlogic = frontVelLoc_interpM>=0.03 & frontVelLoc_interpM<=0.13;
all_CfSmth = conv2(all_Cf,ones(5,1)/5,'same');
CselectLogical = (all_CfSmth>=100 & all_CfSmth<=1250);
all_CfSelect = all_CfSmth;
all_CfSelect(~CselectLogical) = nan;

all_XcSmth = conv2(all_Xc,ones(5,1)/5,'same');
all_XcSelect = all_XcSmth;
all_XcSelect(~CselectLogical) = nan;
figure; plot(all_CfSelect(Xlogic,:),all_XcSelect(Xlogic,:),'.-');




% XcFig = figure; hold on;
% plot([1 1]*Cr,[0 4e-3],'--','Color',[0.2 0.2 0.2],'LineWidth',2);
% text(Cr+10,1.5e-3,'C_R','Color',[0.2 0.2 0.2]);
% plot([1 1]*Cs,[0 4e-3],'--','Color',[0.8 0.8 0.8],'LineWidth',2);
% text(Cs+10,1.5e-3,'C_S','Color',[0.8 0.8 0.8]);
% colormap(Cols4Locs)
% Ch = colorbar;
% Ch.Ticks = movmean(Ch.Ticks,2,'Endpoints','discard');
% TL = cellfun(@num2str,num2cell(findLocs),'UniformOutput',0);
% Ch.TickLabels = TL;
% Ch.Label.String = 'location on block [m]';
