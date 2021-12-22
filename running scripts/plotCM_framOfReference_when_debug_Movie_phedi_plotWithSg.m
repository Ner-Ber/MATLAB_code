plotDamaged=0;
phedis_shiftedRight = bsxfun(@plus,PhediLocation,PhediData.measuredPhedisFromPlot(:)')+shiftByLocation;

center_of_mass = DataStruct.AsperityData.center_of_mass/DataStruct.ExperimentData.res+shiftByLocation;
%--- smooth the cm
N = 4;
center_of_mass = conv2(center_of_mass,ones(N,1)/N,'same');

CM_time = DataStruct.AsperityData.timeCount;

%--- first time of phedi
phedi_t0 = PhediData.timeVec(1);
%--- closest time on CM:
[~,It0] = min(abs(CM_time-phedi_t0));
relevantCM_idxs = It0:(It0+length(PhediData.timeVec)-1);
relevantCM = center_of_mass(relevantCM_idxs,:);
CM_comparison = center_of_mass(It0,:);

figure; hold on;
for i = 1:size(PhediLocation,2)
    if PhediData.slopeIncline(i)>0
        loc_diff = (CM_comparison-phedis_shiftedRight(1,i));
        loc_diff(loc_diff<0) = inf;
        [~,minIdx] = min(loc_diff);
        displacement = phedis_shiftedRight(:,i)-relevantCM(:,minIdx);
        displacement_shift = bsxfun(@minus,displacement,mean(displacement(1:20)));
        plot(PhediData.timeVec,displacement_shift,...
            '*-','Color',FigColors(i,:));
    elseif PhediData.slopeIncline(i)<0
        loc_diff = (CM_comparison-phedis_shiftedRight(1,i));
        loc_diff(loc_diff>0) = inf;
        [~,minIdx] = min(loc_diff);
        displacement = phedis_shiftedRight(:,i)-relevantCM(:,minIdx);
        displacement_shift = bsxfun(@minus,displacement,mean(displacement(1:20)));
        plot(PhediData.timeVec,displacement_shift,...
            '.-','Color',FigColors(i,:));
    else
        if plotDamaged
            plot(phedis_shiftedRight(:,i), PhediData.timeVec,...
                'r--');
        end
    end
end