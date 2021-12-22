%% compare cohesive zone from A(x,t) and from phedis

%% set parameters and stuff
ExpNum = 1;
events = [6];
Colors = MyVaryColor(length(events));
%% create dummies for legend
figure;
hold on;
h = [];
LGND = {};
for EvIdx = 1:length(events)
    h(EvIdx) = plot(nan,nan,'-','Color',Colors(EvIdx,:));
    LGND{EvIdx} = ['event ',num2str(events(EvIdx))];
end
h(EvIdx+1) = plot(nan,nan,'*','Color','k');
LGND{EvIdx+1} = 'from BigPic';
h(EvIdx+2) = plot(nan,nan,'d','Color','k');
LGND{EvIdx+2} = 'from phedis';

%% calc and plot the cohesive
PhediStructCELL = {};
for EvIdx = 1:length(events)
    PhediStructCELL{EvIdx} = Movie_phedi_from_folder_2_data(ExpNum,events(EvIdx));
    cohesiveStruct = IDT_calcCohesiveFromROT(PhediStructCELL{EvIdx}.BigPicRotStruct,[0.06 0.09]);
    plot(cohesiveStruct.Xc,'*','Color',Colors(EvIdx,:));
    
    x_c = phedi_calcCohesiveZoneByVel(PhediStructCELL{EvIdx}.PhediData.PhediVelocity,PhediStructCELL{EvIdx}.PhediData.x_mins_x_tip_4vel);
    plot(x_c,'d','Color',Colors(EvIdx,:));
    title({'X_c comparison',[PhediStructCELL{1}.BigPicRotStruct.details,' event=',num2str(events(EvIdx))]});
end
ylim([-1 1]*1e-2)
legend(h,LGND);