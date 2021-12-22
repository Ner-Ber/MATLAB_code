Event=11;
phedisLoc4Uxx = bsxfun(@plus, analyzePhediCell{Event}.PhediLocation,analyzePhediCell{Event}.measuredPhedisFromPlot');
phedisNums = [3 6 7 10 11]; %[3 6 7 10 11 13 14 15 17 18 19 20 21];
phedisSlopesFiction = zeros(1,size(analyzePhediCell{Event}.PhediLocation,2));
phedisSlopesFiction(phedisNums) = 1;
% phedisSlopesFiction(1:length(phedisSlopesFiction)) = 1;

[UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = phedi_calcUxx(phedisLoc4Uxx,[],[],0.3e-3,phedisSlopesFiction,1);
%% plot spcific phedis and strain
for ii = 1:size(phediPairs,2)
    figure;
    subplot(2,1,1);
    xImagesc = [analyzeAsperityCell{Event}.spatialVec(1) analyzeAsperityCell{Event}.spatialVec(end)];
    yImagesc = [analyzeAsperityCell{Event}.timeCount(1), analyzeAsperityCell{Event}.timeCount(end)];
    imagesc(xImagesc,yImagesc,analyzeAsperityCell{Event}.RowOverTime);
    hold on;
    j1 = phediPairs(1,ii);
    j2 = phediPairs(2,ii);
    plot(phedisLoc4Uxx(:,j1),...
        analyzePhediCell{Event}.timeVec,'r','Linewidth',2);
    plot(phedisLoc4Uxx(:,j2),...
        analyzePhediCell{Event}.timeVec,'g','Linewidth',2);
    hold off;
    
    subplot(2,1,2);
    phedisLoc4UxxZeroed = bsxfun(@minus,phedisLoc4Uxx(:,phediPairs(:,ii)),phedisLoc4Uxx(1,phediPairs(:,ii)));
    yyaxis left
    hold on;
    plot(analyzePhediCell{Event}.timeVec,phedisLoc4UxxZeroed(:,1),'r','Linewidth',2);
    plot(analyzePhediCell{Event}.timeVec,phedisLoc4UxxZeroed(:,2),'g','Linewidth',2);
    yyaxis right
    strechFac = max(abs(UxxFromPhedi(:,ii)))./max(abs(phedisLoc4UxxZeroed(:)));
    plot(analyzePhediCell{Event}.timeVec,UxxFromPhedi(:,ii)/strechFac,'k','Linewidth',1.5);
    title(['UxxFromPhedi multip by ',num2str(1/strechFac)]);
    hold off;
end

% %% plot all strains calculated
% figure;
% hold on;
% for ii = 1:size(phediPairs,2)
%     phedisLoc4UxxZeroed = bsxfun(@minus,phedisLoc4Uxx(:,phediPairs(:,ii)),phedisLoc4Uxx(1,phediPairs(:,ii)));
%     plot(analyzePhediCell{Event}.timeVec,UxxFromPhedi(:,ii),'Linewidth',1.5);
% end
% hold off;