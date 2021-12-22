function chsvFromContact = ROT_findXcFromContact_fromPix(BigPicRotStruct,PhotoLocation)
    
    
    smoothMask = 10;
    [~,I] = min(abs(BigPicRotStruct.x-PhotoLocation));
    rotPhotoCol = BigPicRotStruct.DataMat(:,I);
    rotPhotoColSmth = smooth(rotPhotoCol,smoothMask);
    t_tip = BigPicRotStruct.frontTime_interp(I)/BigPicRotStruct.fps;
    t_mins_t_tipRough= BigPicRotStruct.t - t_tip;
    baseSig = rotPhotoColSmth(t_mins_t_tipRough<-1e-4);
    A0 = mean(baseSig);
    DA0 = 3*var(baseSig);
    %-- fix t_tip
    [Tfix,~] = intersections(t_mins_t_tipRough,rotPhotoColSmth,[-10,10],[1 1]*(A0-DA0));
    TfixMin = min(abs(Tfix));
    t_mins_t_tip = t_mins_t_tipRough+TfixMin;
    
    
%         [x_min_Ttip, relevantLocations] = signal_ChangeTime2Space_mapping(t_mins_t_tip, PhotoLocation, BigPicRotStruct);
%         [x_min_Ttip, relevantLocations] = signal_ChangeTime2Space(t_mins_t_tip, PhotoLocation, BigPicRotStruct);
    
    
    %-- get Cf at xtip
    [~,iLoc] = min(abs(smooth(BigPicRotStruct.frontVelLoc_interpM,5)-PhotoLocation));
    Cf = BigPicRotStruct.frontVel_interpMperS(iLoc);
    x_min_Ttip = -Cf*t_mins_t_tip;
    
    
    bottomLoc = -0.005;
    Yfall = 0.6;
    [~,IxBotm] = min(abs(x_min_Ttip-bottomLoc));
    rotPhotoColShifted = rotPhotoColSmth-rotPhotoColSmth(IxBotm);
    rotPhotoColNorm = rotPhotoColShifted./mean(rotPhotoColShifted(1:50));
    %-- find where the smoothed curve meets 1:
    [XDelta,~] = intersections(x_min_Ttip,rotPhotoColNorm,[-100,100],[1 1]*0.95);
    %     DX = min(XDelta(XDelta>0));
    [~,DXi] = min(XDelta);
%     x_min_TtipFix = x_min_Ttip-XDelta(DXi);
    x_min_TtipFix = x_min_Ttip;
    
    %-- find where the falls is of wanted magnitude:
    [XcPosble,~] = intersections(x_min_TtipFix,rotPhotoColNorm,[-100,100],[1 1]*Yfall);
    [~,Xci] = min(abs(XcPosble));
    Xc = XcPosble(Xci);
    [TcPosble,~] = intersections(t_mins_t_tip,rotPhotoColNorm,[-100,100],[1 1]*Yfall);
    [~,ti] = min(abs(TcPosble-t_tip));
    Tc = TcPosble(ti);
    DA = mean(rotPhotoColShifted(1:50));
    [~,IxBotmFix] = min(abs(x_min_TtipFix-bottomLoc));
    A0 = DA+rotPhotoColSmth(IxBotmFix);
    
    
    if isempty(Xc)
        Xc = nan;
    end
    chsvFromContact.Xc = Xc;
    chsvFromContact.loc = PhotoLocation;
    chsvFromContact.Tc = Tc;
    chsvFromContact.DA = DA;
    chsvFromContact.A0 = A0;
    chsvFromContact.Cf = Cf;
    chsvFromContact.Dx = XDelta(DXi);
    chsvFromContact.Dt = TfixMin;
    chsvFromContact.Yfall = Yfall;
    chsvFromContact.bottomLoc = bottomLoc;
    chsvFromContact.x_min_TtipFix = x_min_TtipFix;
    chsvFromContact.t_mins_t_tip = t_mins_t_tip;
    chsvFromContact.rotPhotoColNorm = rotPhotoColNorm;
    
    
    %%
    %     [~,Icf] = min(abs(BigPicRotStruct.frontVelLoc_interpM-PhotoLocation));
    %     Cfloc = BigPicRotStruct.frontVel_interpMperS(Icf);
    %     figure(XcFig);
    %     P = plot(x_min_TtipFix,rotPhotoColNorm,'.-','DisplayName',[BigPicRotStruct.details, 'Cf=',num2str(Cfloc)]);
    %     plot(Xc,Yfall,'o','LineWidth',2,'Color',P.Color);
    
    
end