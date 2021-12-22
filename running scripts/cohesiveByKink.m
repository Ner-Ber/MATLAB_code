%% calc cohesive length
% DirNames = dir_names;
FileNames = {'allPhediStructures_exper_1.mat','allPhediStructures_exper_2.mat', 'allPhediStructures_exper_3.mat'};
X_kink_Cell= cell(size(FileNames));
cf_Cell = cell(size(FileNames));
Ev_Cell = cell(size(FileNames));
for d = 1:length(FileNames)
    %--- load the data
    S = load(FileNames{d});
    PhediStruct_d = S.PhediStruct;
    %--- read and analyze data of events
    X_kink_VEC = [];
    cf_VEC = [];
    Ev_VEC = [];
    for Ev = 1:length(PhediStruct_d)
        if isempty(PhediStruct_d{Ev})
            continue
        end
        [X_kink,~] = phedi_findCohesiveZoneKink(PhediStruct_d{Ev}.PhediData);
        %         X_kink = abs(mean(X_kink));
        X_kink_VEC = cat(1,X_kink_VEC,X_kink);
        cf_VEC = cat(1,cf_VEC,repmat(PhediStruct_d{Ev}.BigPicRotStruct.AvgCf_scratches,length(X_kink),1));
        Ev_VEC = cat(1,Ev_VEC,repmat(Ev,length(X_kink),1));
    end
    X_kink_Cell{d} = X_kink_VEC;
    cf_Cell{d} = cf_VEC;
    Ev_Cell{d} = Ev_VEC;
end
%% plot
markerList = {'o','+','*','x','s','d','^','v','>','<','p','h'};
colorList = {rgb('Blue'),rgb('Red'),rgb('Orange'),rgb('Purple'),rgb('Green'),...
    rgb('HotPink'), rgb('Yellow'),rgb('Gray')};
figure; hold on;
for d = 1:length(FileNames)
    U_ev = unique(Ev_Cell{d});
    for u_i = U_ev(:)'
        Ev_logical = Ev_Cell{d}==u_i;
        Mean_X_kink = mean(abs(X_kink_Cell{d}(Ev_logical)));
        Min_X_kink = min(abs(X_kink_Cell{d}(Ev_logical)));
        Max_X_kink = max(abs(X_kink_Cell{d}(Ev_logical)));
        cf_this = unique(cf_Cell{d}(Ev_logical));
%         errorbar(cf_this,Mean_X_kink,Min_X_kink,Max_X_kink,'*','Color',colorList{d},'LineWidth',1.2);
        plot(cf_this,Mean_X_kink,'o','Color',colorList{d},'LineWidth',1.2,'MarkerFaceColor',colorList{d});
    end
    
end


