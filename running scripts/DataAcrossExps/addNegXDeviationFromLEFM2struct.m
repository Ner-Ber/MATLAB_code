D = dir_names
k1 = strfind(D,'.mat')
k2 = strfind(D,'allPhedi')
F = find(~cellfun(@isempty,k1) & ~cellfun(@isempty,k2))
A_avgDist = 0.01;
for j=1:length(F)
    %% iterate on events
    try
        disp(['loading ',D{F(j)}]);
        load(D{F(j)});
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
        relevantEvents = relevantEvents(relevantEvents~=1);
        
        [maxXi_loc_vec,maxYi_loc_vec,maxXi_vel_vec,maxYi_vel_vec] = ...
            deal(zeros(length(relevantEvents),1));
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                
                [maxXi_loc,maxYi_loc,maxXi_vel,maxYi_vel] =...
                    phedi_intersections_w_LEFM(PhediStructCell{relevantEvents(iEv)},'00');
                
                [maxXi_loc_vec(iEv),maxYi_loc_vec(iEv),maxXi_vel_vec(iEv),maxYi_vel_vec(iEv)] = ...
                    deal(maxXi_loc,maxYi_loc,maxXi_vel,maxYi_vel);
                
            catch
            end
        end
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        J = arrayfun(@(S) strcmpi(S.Name,Name), myStruct);
        
        myStruct(J).maxXi_loc_vec_allSlps = maxXi_loc_vec;
        myStruct(J).maxYi_loc_vec_allSlps = maxYi_loc_vec;
        myStruct(J).maxXi_vel_vec_allSlps = maxXi_vel_vec;
        myStruct(J).maxYi_vel_vec_allSlps = maxYi_vel_vec;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end