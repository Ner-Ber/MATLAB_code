function Cf = sg_calc_velocity_from_sg(sg_data_struct, relevant_SGs)
% sg_calc_velocity_from_sg returns the front speed as averages between all
%combinations of SGs in the input. Uses the minimal point of the Uxx
%component as an indicator of fron location.
%Cf is returns as an absolute value in m/s. 


%-- -use more than 1 sg:
if length(relevant_SGs)<2
    if relevant_SGs==1
        relevant_SGs = [relevant_SGs relevant_SGs+1];
    else
        relevant_SGs = [relevant_SGs relevant_SGs-1];
    end
end
numSG_2_calc = length(relevant_SGs);

relevantUxx = sg_data_struct.Uxx(:,relevant_SGs);
relevantSgLocation = sg_data_struct.x_sg(relevant_SGs);

[~,IDXs] = min(relevantUxx,[],1);



AllPerms = nchoosek(1:numSG_2_calc,2);
allCalcCFs = zeros(size(AllPerms,1),1);
for permIdx = 1:size(AllPerms,1)
    TIMES = sg_data_struct.t(IDXs(AllPerms(permIdx,:)));
    Locs = relevantSgLocation(AllPerms(permIdx,:));
    allCalcCFs(permIdx) = abs(diff(Locs)/diff(TIMES))*10;
end
Cf = mean(allCalcCFs);

end