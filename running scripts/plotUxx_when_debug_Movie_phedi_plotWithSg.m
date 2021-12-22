phedisLoc4Uxx = bsxfun(@plus, PhediData.PhediLocation,PhediData.measuredPhedisFromPlot');
phedisSlopesFiction = zeros(1,size(PhediData.PhediLocation,2));
phedisSlopesFiction(1:length(phedisSlopesFiction)) = 1;

[UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx(phedisLoc4Uxx,[],[],0.15e-3,phedisSlopesFiction,1);

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