zeroThresh = 5e-8;
relevantEvents = find(~cellfun(@isempty,PhediStruct));
% relevantEvents = EventsWGoodPhds;
AllCohsvLength = cell(size(relevantEvents));
doPlot = 0;
Cfs = zeros(size(PhediStruct));
Gamma = 11;
for ev = 6 % 1:length(relevantEvents)
    DataStructUpdated = Movie_fix_t_tips(PhediStruct{relevantEvents(ev)});
    relevantPhed = find(DataStructUpdated.PhediData.slopeIncline==1);
    AllCohsvLength{ev} = zeros(1,max(relevantPhed));
    Cfs(ev) = DataStructUpdated.PhediData.Cf;
    PhediLocReduced = phedi_reduceSinFromLocation(DataStructUpdated.PhediData);
    t_mins_t_tip = bsxfun(@minus,DataStructUpdated.PhediData.timeVec,DataStructUpdated.PhediData.t_tips(:)');
    x_min_x_tip = -Cfs(ev)*t_mins_t_tip;
    sol = CrackSolutionForh_GammaChange(Cfs(ev)/1255,-0.5,1e-9,Gamma,-0.2:1e-7:0.8);
    [power,sqrtPref] = Crack_UxFromInterfaceSolution(sol);
    x = sol.x;
    y = sqrtPref*(x.^power);
    y((-x)>0) = 0;
    [x,Order] = sort(-x);
    y = y(Order);
    
    figure;
    hold on;
    plot(commonX,y_c,'k','LineWidth',1.2);
    PhediColors = MyVaryColor(length(DataStructUpdated.PhediData.slopeIncline));
    LGND = {'LEFM fit'};
    
    for phd = relevantPhed
        [~, shifted_xAxis] = phedi_findShiftByFit2PowerLaw(x_min_x_tip(:,phd),PhediLocReduced(:,phd),power,sqrtPref, [-0.025 0.015]);
        [commonX, PhediLoc_c, y_c] = commonVectorPoints(shifted_xAxis,PhediLocReduced(:,phd),x,y,1e-5);
        modelDiff = PhediLoc_c - y_c;
        [x0,y0] = intersections([-100 100],[0 0],commonX,modelDiff);
        [closestx0,Ix0min] = min(abs(x0));
%         closestx0 = closestx0*sign(x0(Ix0min));
%         farCohesiveEdge = x(find(commonX>0 & abs(modelDiff)<zeroThresh, 1, 'first'));
%         farCohesiveEdge = x0(Ix0min+1);
%         cohesiveZone = farCohesiveEdge - closestx0;
%         AllCohsvLength{ev}(phd) = cohesiveZone;
        
        hold on; plot(commonX-0.004, PhediLoc_c,'*-','Color',PhediColors(phd,:),'MarkerSize',0.5);
%         plot(commonX(commonX>closestx0 & commonX<farCohesiveEdge),PhediLoc_c(commonX>closestx0 & commonX<farCohesiveEdge),'r','LineWidth',2);
        title(['Ev=',num2str(ev),'    C_f=',num2str(Cfs(ev))]);
        LGND = cat(1,LGND,['phd=',num2str(phd)]);
    end
    legend(LGND);
end

% figure; hold on;
% for i=1:length(relevantEvents)
% plot(Cfs_mod(i),mean(AllCohsvLength{i}),'o','MarkerFaceColor','m');
% end





% % relevantEvents = find(~cellfun(@isempty,PhediStruct));
% relevantEvents = EventsWGoodPhds;
% AllCohsvLength = cell(size(relevantEvents));
% doPlot = 0;
% Cfs = zeros(size(PhediStruct));
% for ev = 1:length(relevantEvents)
%     DataStructUpdated = Movie_fix_t_tips(PhediStruct{relevantEvents(ev)});
%     Movie_phedi_plotWithSg(DataStructUpdated,1,'sep',{'loc'});
%     relevantPhed = find(DataStructUpdated.PhediData.slopeIncline==1);
%     AllCohsvLength{ev} = zeros(1,max(goodPhedis{ev}));
%     Cfs(ev) = DataStructUpdated.PhediData.Cf;
%     for phd = goodPhedis{ev}
%         [cohesiveZone, T_mean,PhediLoc_mean, newPowerStruct] = phedi_CohesiveFromFit(DataStructUpdated,[],[],[],[],phd,doPlot);
%         if doPlot
%             title([DataStructUpdated.BigPicRotStruct.details,' phedi=',num2str(phd)]);
%         end
%         AllCohsvLength{ev}(phd) = abs(diff(cohesiveZone));
%     end
% end
% 
% figure; hold on;
% for i=1:length(relevantEvents)
% plot(Cfs_mod(i),mean(AllCohsvLength{i}),'o','MarkerFaceColor','m');
% end