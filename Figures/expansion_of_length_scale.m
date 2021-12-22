%% plot expansion of length scale
[Cd, Cs, Cr, nu, ro, E, mu, ~, PlaneStrain, tau_p, Xc0] = CrackSolutionMaterialProperties;
% A1 = load('G:\Frics\2018-10-10\DAandREsidualsStruct.mat');
% A2 = load('G:\Frics\2018-8-29\DAandREsidualsStruct.mat');
A1 = load('E:\Frics\2018-10-10\DAandREsidualsStruct.mat');
A2 = load('E:\Frics\2018-8-29\DAandREsidualsStruct.mat');
myStruct20181010 = A1.myStruct;
myStruct20180829 = A2.myStruct;
myStruct = myStruct20180829;
myStruct((length(myStruct20180829)+1):(length(myStruct20180829)+length(myStruct20181010))) = myStruct20181010(:);
%--- clear empty fields
emptyLogic = arrayfun(@(A) isempty(A.Cf_vec) , myStruct);
myStruct(emptyLogic)=[];

%--- set velocity scale
Cf_scale = linspace(0,1260,9);
CfVals = movmean(Cf_scale,2,'Endpoints','discard');
Colors = MyVaryColor(length(CfVals),flipud(My_colorMap));

%--- create cell of event details
Events_vec = arrayfun(@(S) S.Events,myStruct,'UniformOutput',0);
Events_vec = cell2mat(Events_vec);
Events_vec = Events_vec(:);
Events_cell = num2cell(Events_vec);
Names_cell = arrayfun(@(S) repmat({S.Name},length(S.Events),1),myStruct,'UniformOutput',0);
Names_cell = Names_cell(:);
Names_cell_combined = {};
for i=1:length(Names_cell)
    Names_cell_combined = cat(1,Names_cell_combined,Names_cell{i});
end
event_details = cellfun(@(A,B) [A,' ',num2str(B)],Names_cell_combined,Events_cell,'UniformOutput',0);

%--- define damaged events
damaged_events20180829 = {...           % for 29-8-2018
    '12-8-52 15','12-8-52 17',...
    '14-2-11 6','14-2-11 12',...
    '14-45-14 10','14-45-14 19','14-45-14 28','14-45-14 12',...
    '15-9-25 19',...
    '15-28-23 29','15-28-23 19','15-28-23 18','15-28-23 23','15-28-23 13','15-28-23 8',...
    '18-1-9 14','18-1-9 11',...
    '18-20-0 24','18-20-0 19','18-20-0 7','18-20-0 4'...
    };
damaged_events_wDate20180829 = cellfun(@(A) ['2018-8-29 ',A],damaged_events20180829,'UniformOutput',0);

damaged_events20181010 = {...           % for 10-10-2018
    '12-2-22 28','12-2-22 20',...
    '14-28-36 4',...
    '14-48-1 12','14-48-1 21',...
    '20-28-52 14',...
    '20-48-22 16','20-48-22 3','20-48-22 5',...
    };
damaged_events_wDate20181010 = cellfun(@(A) ['2018-10-10 ',A],damaged_events20181010,'UniformOutput',0);

damaged_events = damaged_events_wDate20180829(:);
damaged_events(length(damaged_events_wDate20180829)+(1:length(damaged_events_wDate20181010))) = damaged_events_wDate20181010(:);
damaged_events_logic = ~ismember(event_details,damaged_events);

N = cellfun(@(S) min(strfind(S,' ')), event_details,'UniformOutput',0);
event_dates = cellfun(@(A,n) A(1:(n-1)),event_details,N,'UniformOutput',0);
[~,~,groupCatagrs0] = unique(event_dates);
presentedData = cellfun(@(A) ~isempty(strfind(A,'17-18-6')),Names_cell_combined);
groupCatagrs = groupCatagrs0;
groupCatagrs(presentedData) = 3;
%--- extract Cf
Cf_Cell = arrayfun(@(S) S.Cf_vec,myStruct,'UniformOutput',0);
Cf_Cell = Cf_Cell(:);
Cf_mat = cell2mat(Cf_Cell);
Cf_logical = Cf_mat<=1260 & Cf_mat>=0;

%--- define display logical
Logical_tot = Cf_logical & damaged_events_logic;

Cf_filter = Cf_mat(Logical_tot);

%--- exctract peak velocity
Vpeak_All_cell = arrayfun(@(S) S.peakVel_vecAll,myStruct,'UniformOutput',0);
Vpeak_All_cell = Vpeak_All_cell(:);
Vpeak_All = cell2mat(Vpeak_All_cell);
%--- extract gamma
G_Cell = arrayfun(@(S) S.Gamma,myStruct,'UniformOutput',0);
G_Cell = G_Cell(:);
G_mat = cell2mat(G_Cell);
G_filter = G_mat(Logical_tot);

%--- tau peak by model
Vpeak_All_filter = Vpeak_All(Logical_tot);
tau_p_vec = nan(size(G_filter));
peakVel_vec = nan(size(G_filter));
for i=1:length(G_filter)
    if Cf_filter(i)>=Cr
        thisCf = Cr-1;
    else
        thisCf = Cf_filter(i);
    end
    [X_peak, peakVel_vec(i)] = cohesive_findPeakOfVelByIterations(thisCf/Cr,1e-9,G_filter(i),'L',L*10,'exp');
    tau_p_vec(i) = CohesiveSolutionGetTau_p(thisCf/Cr,G_filter(i),L*10,'exp');
end
Ratio = tau_p_vec./peakVel_vec;
tau_p_semiPredictor = Vpeak_All_filter.*Ratio;
slmRatVsCf = cohesive_tau_p_predictorByVx('exp');
Ratio_4Cf = ppval(Cf_filter,slmRatVsCf);
tau_p_byModel = Ratio_4Cf.*Vpeak_All_filter;


%--- prepare and plot
Vpeak_All_filter = Vpeak_All(Logical_tot);
groupCatagrs_filter = groupCatagrs(Logical_tot);

%% calculate Xc:

Xc = 2.*G_filter.*Cf_filter./tau_p_byModel./Vpeak_All_filter;

figure; plot(Cf_filter/Cr,Xc,'o');
pbaspect([4,3,1])
