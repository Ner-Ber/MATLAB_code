%% create figure of system

%% load for sasmples
% load('G:\Frics\2018-10-10\allPhediStructures_165017.mat');
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
for iEv = 1:length(relevantEvents)
    PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
    PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
end

EvNum = 21;
FontSize = 10;


%% create figure and axes
F1 = figure;
F1.Color=[1 1 1];
F1.WindowStyle = 'normal';
F1.Units = 'centimeters';
F1.Position([3,4]) = [18 8];
systemDescprition = axes;
systemDescprition.Position = [0.13, 0.6, 0.8, 0.35];   % top
smallPic_ROT = axes;
smallPic_ROT.Position = [0.13, 0.12, 0.35, 0.35];   % bottom left
bigPic_ROT = axes;
bigPic_ROT.Position = [0.58, 0.12, 0.35, 0.35];   % bottom right

%% system description
systemDescprition.XAxis.TickLabels=[];
systemDescprition.YAxis.TickLabels=[];
systemDescprition.XAxis.TickValues=[];
systemDescprition.YAxis.TickValues=[];
%% BigPic Row over time
BigPicROT_fig = bigPic_ROT;
IDT_PlotRowOverTime(PhediStructCell{EvNum}.BigPicRotStruct);
axis([[-0.001 0.151]*1 [-2 18]*1e-4]);
caxis([0.7 1.03]);
for i=1:length(BigPicROT_fig.Children)
    if  ~strcmpi(BigPicROT_fig.Children(i).Type,'image')
        BigPicROT_fig.Children(i).Visible='off';
    end
end

bigPic_ROT.Title.Visible='off';
bigPic_ROT.YAxis.TickValues = [0:0.5:1.5]*1e-3;
bigPic_ROT.YAxis.TickLabels = [0:0.5:1.5];
ylabel('t [ms]');

%-- renovate colorbar

F1.Children(1).Label.String = 'A/A_0';
F1.Children(1).Position(3) = 0.018;
F1.Children(1).Position(1) = 0.9;
F1.Children(1).Label.Position(1) = 2.5;
bigPic_ROT.Position = [0.58, 0.12, 0.3, 0.35];


%% small pic ROT
lineWidth = 3;

SmallPicROT_fig = smallPic_ROT;
Movie_phedi_plotWithSg(PhediStructCell{EvNum},'00',{'ROT'},[],[],[],SmallPicROT_fig);
axis([[0.07875 0.0807]*1 [-1.5 4]*1e-3]);
HideLines = [arrayfun(@(S) strcmpi(S.Type,'line') & isempty(S.DisplayName), SmallPicROT_fig.Children(1:end-1))];
for i=find(HideLines(:)'); SmallPicROT_fig.Children(i).Visible='off'; end
for i=find(~HideLines(:)')
    smallPic_ROT.Children(i).Marker = 'none';
    smallPic_ROT.Children(i).LineWidth = 3;
    smallPic_ROT.Children(i).Color = rgb('DeepSkyBlue');
end
colormap(smallPic_ROT,'gray');
ylabel('t [ms]');
xlabel('x [m]');
smallPic_ROT.Title.Visible='off';
smallPic_ROT.YAxis.TickValues = [-1:3]*1e-3;
smallPic_ROT.YAxis.TickLabels = [-1:3];



%% set all sub figure
set(findall(F1,'-property','FontSize'),'FontSize',FontSize);