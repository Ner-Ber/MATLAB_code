function chsvFromContact = ROT_findXcFromContact_fromFrame(BigPicRotStruct,PhotoLocation)
    
    
%     smoothMask = 10;
%     [~,I] = min(abs(BigPicRotStruct.x-PhotoLocation));
%     rotPhotoCol = BigPicRotStruct.DataMat(:,I);
%     rotPhotoColSmth = smooth(rotPhotoCol,smoothMask);
%     t_tip = BigPicRotStruct.frontTime_interp(I)/BigPicRotStruct.fps;
%     t_mins_t_tipRough= BigPicRotStruct.t - t_tip;
%     baseSig = rotPhotoColSmth(t_mins_t_tipRough<-1e-4);
%     A0 = mean(baseSig);
%     DA0 = 3*var(baseSig);
%     %-- fix t_tip
%     [Tfix,~] = intersections(t_mins_t_tipRough,rotPhotoColSmth,[-10,10],[1 1]*(A0-DA0));
%     TfixMin = min(abs(Tfix));
%     t_mins_t_tip = t_mins_t_tipRough+TfixMin;
%     
%     
%     %         [x_min_Ttip, relevantLocations] = signal_ChangeTime2Space_mapping(t_mins_t_tip, PhotoLocation, BigPicRotStruct);
%     %         [x_min_Ttip, relevantLocations] = signal_ChangeTime2Space(t_mins_t_tip, PhotoLocation, BigPicRotStruct);
%     
%     
%     %-- get Cf at xtip
%     [~,iLoc] = min(abs(smooth(BigPicRotStruct.frontVelLoc_interpM,5)-PhotoLocation));
%     Cf = BigPicRotStruct.frontVel_interpMperS(iLoc);
%     x_min_Ttip = -Cf*t_mins_t_tip;
%     
%     
%     bottomLoc = -0.01;
%     Yfall = 0.6;
%     [~,IxBotm] = min(abs(x_min_Ttip-bottomLoc));
%     rotPhotoColShifted = rotPhotoColSmth-rotPhotoColSmth(IxBotm);
%     rotPhotoColNorm = rotPhotoColShifted./mean(rotPhotoColShifted(1:50));
%     %-- find where the smoothed curve meets 1:
%     [XDelta,~] = intersections(x_min_Ttip,rotPhotoColNorm,[-100,100],[1 1]*0.95);
%     %     DX = min(XDelta(XDelta>0));
%     [~,DXi] = min(XDelta);
%     %     x_min_TtipFix = x_min_Ttip-XDelta(DXi);
%     x_min_TtipFix = x_min_Ttip;
%     
%     %-- find where the falls is of wanted magnitude:
%     [XcPosble,~] = intersections(x_min_TtipFix,rotPhotoColNorm,[-100,100],[1 1]*Yfall);
%     [~,Xci] = min(abs(XcPosble));
%     Xc = XcPosble(Xci);
%     [TcPosble,~] = intersections(t_mins_t_tip,rotPhotoColNorm,[-100,100],[1 1]*Yfall);
%     [~,ti] = min(abs(TcPosble-t_tip));
%     Tc = TcPosble(ti);
%     DA = mean(rotPhotoColShifted(1:50));
%     [~,IxBotmFix] = min(abs(x_min_TtipFix-bottomLoc));
%     A0 = DA+rotPhotoColSmth(IxBotmFix);
%     
%     
%     if isempty(Xc)
%         Xc = nan;
%     end
%     chsvFromContact.Xc = Xc;
%     chsvFromContact.Tc = Tc;
%     chsvFromContact.DA = DA;
%     chsvFromContact.A0 = A0;
%     chsvFromContact.Cf = Cf;
%     chsvFromContact.Dx = XDelta(DXi);
%     chsvFromContact.Dt = TfixMin;
%     chsvFromContact.Yfall = Yfall;
%     chsvFromContact.bottomLoc = bottomLoc;
%     chsvFromContact.x_min_TtipFix = x_min_TtipFix;
%     chsvFromContact.t_mins_t_tip = t_mins_t_tip;
%     chsvFromContact.rotPhotoColNorm = rotPhotoColNorm;
    
    %% by single frame
    %-- find frame
    smoothMask = 10;
    Yfall = 0.6;
    frontTime = BigPicRotStruct.frontTime_interp/BigPicRotStruct.fps;       % time vector for front
    [X0,T0] = intersections(BigPicRotStruct.x,frontTime,PhotoLocation*[1 1],[-1,1]);    % get time when from reches 'PhotoLocation'
    [~,I] = min(abs(BigPicRotStruct.t-T0));
    X_minXtip = BigPicRotStruct.x-X0;
    snapsht = BigPicRotStruct.DataMatNorm(I,:);   %snapshot relevant to requested location 'PhotoLocation'
    bottomLoc = -0.005;
    [~,BottmPointPix] = min(abs(X_minXtip-bottomLoc));
    deltaM = 1e-3; % define the width of window on which to measure the bottom value
    deltaP = round(deltaM/BigPicRotStruct.res);
    BottmPointVal = mean(snapsht(BottmPointPix-deltaP:BottmPointPix+deltaP));
    snapshtStrchd = (snapsht-BottmPointVal)./(1-BottmPointVal);
    snapshtSmth = smooth(snapshtStrchd,smoothMask);
    [XcPosble,~] = intersections(X_minXtip,snapshtSmth,[-100,100],Yfall*[1 1]);
    Xc = max(XcPosble(XcPosble<0));
    if isempty(Xc)
        Xc = nan;
    end

    
    [~,iLoc] = min(abs(smooth(BigPicRotStruct.frontVelLoc_interpM,5)-PhotoLocation));
    Cf = BigPicRotStruct.frontVel_interpMperS(iLoc);
    
    Tc = -Xc./Cf;
    
    [~,zoerPointPix] = min(abs(X_minXtip+bottomLoc));
    A0 = mean(BigPicRotStruct.DataMat(I,zoerPointPix-deltaP:zoerPointPix+deltaP));
    DA = A0 - mean(BigPicRotStruct.DataMat(I,BottmPointPix-deltaP:BottmPointPix+deltaP));
    
    chsvFromContact.Xc = Xc;
    chsvFromContact.loc = PhotoLocation;
    chsvFromContact.Tc = Tc;
    chsvFromContact.DA = DA;
    chsvFromContact.A0 = A0;
    chsvFromContact.Cf = Cf;
    chsvFromContact.Yfall = Yfall;
    chsvFromContact.bottomLoc = bottomLoc;
    chsvFromContact.rotPhotoRowNorm = snapshtSmth;
    
    %%
    %     [~,Icf] = min(abs(BigPicRotStruct.frontVelLoc_interpM-PhotoLocation));
    %     Cfloc = BigPicRotStruct.frontVel_interpMperS(Icf);
    %     figure(XcFig);
    %     P = plot(x_min_TtipFix,rotPhotoColNorm,'.-','DisplayName',[BigPicRotStruct.details, 'Cf=',num2str(Cfloc)]);
    %     plot(Xc,Yfall,'o','LineWidth',2,'Color',P.Color);
    
    
end