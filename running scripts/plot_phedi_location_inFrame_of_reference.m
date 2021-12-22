%% create common time vector
t_mins_t_tip = bsxfun(@minus,PhediStructEv8.PhediData.timeVec,PhediStructEv8.PhediData.t_tips(:)');
t_mins_t_tip_sortd = sort(t_mins_t_tip(:));
%% create location vectors to fit general time vector
PhediLocation = PhediStructEv8.PhediData.PhediLocation;
PhediLocsInterped = nan(length(t_mins_t_tip_sortd),16);
for i=1:16; PhediLocsInterped(:,i)=interp1(t_mins_t_tip(:,i),PhediLocation(:,i),t_mins_t_tip_sortd,'linear','extrap'); end
PhediLocsInterped(PhediLocsInterped==0)=nan;
%% creat mean trajectory
PhediLocsInterpedZeroed = bsxfun(@minus,PhediLocsInterped,mean(PhediLocsInterped(300:400,:),1,'omitnan'));
PhediLocsInterpedMeaned = mean(PhediLocsInterpedZeroed,2,'omitnan');
PhediLocsInterpedStd = std(PhediLocsInterpedZeroed,0,2,'omitnan');
PhediLocsInterpedZeroedSummed = sum(PhediLocsInterpedZeroed,2,'omitnan');
% PhediLocsInterpedStd_presmooth = std(conv2(PhediLocsInterpedZeroed,ones(1,100),100),0,2,'omitnan');
% PhediLocsInterpedP2P = peak2peak(PhediLocsInterpedZeroed,2,'omitnan');
%% plot 
figure;
plot(t_mins_t_tip_sortd,PhediLocsInterpedZeroedSummed);
[~,II] = min(PhediLocsInterpedZeroedSummed);
hold on; plot([1 1]*t_mins_t_tip_sortd(II),[-4 10]*1e-5,'r');
%% create trajectories in center of mass frame of reference
PhedisLocsCM = nan(size(PhediLocation));
for i=1:16; PhedisLocsCM(:,i)= PhediLocation(:,i)-PhediLocsInterpedMeaned(ismember(t_mins_t_tip_sortd,t_mins_t_tip(:,i))); end

%% plot this
figure; plot(t_mins_t_tip,PhedisLocsCM);
%% copy colors and textures from other plots
GCA = get(gca,'Children');
thisFig = get(gca,'Children');
for i=length(thisFig):-1:1
thisFig(i).Color = GCA(i).Color;
thisFig(i).Marker = GCA(i).Marker;
thisFig(i).DisplayName = GCA(i).DisplayName;
end

%% calculate strain in CM FOR
%--- omit nans:
NoNansLogical = ~isnan(sum(PhediLocsInterped,2));
PhediLocsNanOmitd = PhediLocsInterped(NoNansLogical,:);
t_mins_t_tip_NoNans = t_mins_t_tip_sortd(NoNansLogical);
[UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx_ofPairs(PhediLocsNanOmitd);


