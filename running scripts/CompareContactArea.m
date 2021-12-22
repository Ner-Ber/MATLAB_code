%% compare locations of contact area
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
Location = 0.02:0.02:0.14;
Markers = {'*','o','s','d','^','v','x','p','h','+','>','<'};
Colors = distinguishable_colors(23);
% ContactCmp = figure;
figure(ContactCmp);
ax = gca;
hold on;
%-- iterate on events:
for iEv = 1:length(relevantEvents)
    %-- iterate on locations:
    for iLoc = 1:length(Location)        
        BigPicRotStruct = PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct;
        thisLoc = Location(iLoc);
        x_minXtip_4spatial = PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct.x - thisLoc;
        frontTime = BigPicRotStruct.frontTime_interp./BigPicRotStruct.fps;
        [~,I1] = min(abs(x_minXtip_4spatial));
        t_minTtip_4tmprl = BigPicRotStruct.t-frontTime(I1);
        [~,I2] = min(abs(t_minTtip_4tmprl));
        A_spatial = BigPicRotStruct.DataMatNorm(I2,:);
        
        Name = [BigPicRotStruct.details,' loc=',num2str(thisLoc)];
        plot(x_minXtip_4spatial,A_spatial,'Color',Colors(iEv,:),'Marker',Markers{iLoc},'DisplayName',Name)
    end
    %-- restart color order
    ax.ColorOrderIndex = 1;
end
