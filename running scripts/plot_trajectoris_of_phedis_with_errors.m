samples2Use = [3 6 7 10 11 13 14 15 17 18 19 20 21];
LeftSamples = [3 7 11 13 15 17  19  21];
RightSamples = setdiff(samples2Use,LeftSamples);

Event =  11;
PhediLocation = analyzePhediCell{Event}.PhediLocation;
Phedis2Zero = bsxfun(@minus,PhediLocation,PhediLocation(1,:));
Phedis2Zero = bsxfun(@minus,Phedis2Zero,mean(Phedis2Zero(:,samples2Use),2));    % optional row, this will plot the trajectories in the center of mass frame of reference
t = analyzePhediCell{Event}.timeVec;

StatStruct = timeSamplesGetStatistics(Phedis2Zero, samples2Use);
m_tot = StatStruct.mean;
dm_tot = StatStruct.std;

StatStruct = timeSamplesGetStatistics(Phedis2Zero, LeftSamples);
m_left = StatStruct.mean;
dm_left = StatStruct.std;

StatStruct = timeSamplesGetStatistics(Phedis2Zero, RightSamples);
m_right = StatStruct.mean;
dm_right = StatStruct.std;

figure;
hold on;
h = [];
%--- plot errors
h(1) = fill([t;flipud(t)],[m_tot-dm_tot;flipud(m_tot+dm_tot)],[.4 1 1],'linestyle','none','FaceAlpha',0.5);
h(2) = fill([t;flipud(t)],[m_left-dm_left;flipud(m_left+dm_left)],[1 0 0],'linestyle','none','FaceAlpha',0.5);
h(3) = fill([t;flipud(t)],[m_right-dm_right;flipud(m_right+dm_right)],[.2 1 .6],'linestyle','none','FaceAlpha',0.5);

%--- plot means
h(4) = plot(analyzePhediCell{Event}.timeVec,m_tot,'b','LineWidth',1.8);
h(5) = plot(analyzePhediCell{Event}.timeVec,m_left,'r','LineWidth',1.8);
h(6) = plot(analyzePhediCell{Event}.timeVec,m_right,'g','LineWidth',1.8);


legend(h(4:end),{'all','left','right'})