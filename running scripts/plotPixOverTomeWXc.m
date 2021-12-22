figure; hold on;

relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
Cf_re_vec = linspace(0,1255,20);
Cf_cols = MyVaryColor(20);


colormap(Cf_cols)
Ch = colorbar;
Ch.Ticks = movmean(0:20,2,'Endpoints','discard');
TL = cellfun(@num2str,num2cell(Cf_re_vec),'UniformOutput',0);
Ch.TickLabels = TL;
Ch.Label.String = 'Cf@loc [m/s]';

for I=1:length(relevantEvents)
    Cf_vec = nan(length(findLocs),1);
    Xc_vec = nan(length(findLocs),1);
    for f = 1:length(findLocs)
        try
            ChsvFromCntct = ROT_findXcFromContact(PhediStructCell{relevantEvents(I)}.BigPicRotStruct,findLocs(f));
            Xc_i = ChsvFromCntct.Xc;
            Cf_i = ChsvFromCntct.Cf;
            [~,ci] = min(abs(Cf_re_vec-Cf_i));
            Name = [PhediStructCell{relevantEvents(I)}.BigPicRotStruct.details,' loc=',num2str(findLocs(f))];
%             plot(ChsvFromCntct.x_min_TtipFix,ChsvFromCntct.rotPhotoColNorm,'.-','Color',Cf_cols(ci,:),'DisplayName',Name,'LineWidth',1.5);
            plot(ChsvFromCntct.x_min_TtipFix,ChsvFromCntct.rotPhotoColNorm,'.-','Color',Cf_cols(ci,:),'DisplayName',Name,'LineWidth',1.5);
            plot(Xc_i,ChsvFromCntct.Yfall,'o','LineWidth',2,'Color',Cf_cols(ci,:));
        catch
        end
    end
    
end