%chsvFromContact = ROT_findXcFromContact(BigPicRotStruct,varargin)
%[Afall, bottomDist, dx, smt] = setDefaults4function(varargin,0.67,-0.01, 0.0005, 5);


function chsvFromContact = ROT_findXcFromContact(BigPicRotStruct,varargin)
    
    %% set defaults and parameters
    [Afall, bottomDist, dx, smt] = setDefaults4function(varargin,0.67,-0.01, 0.0005, 5);
%     Afall = 0.67;
%     bottomDist = -0.01;
%     dx = 0.005;
    
    %%
    XX = 1:length(BigPicRotStruct.x);
    %-- front's time:
    frontTimePix = BigPicRotStruct.frontTime_interp;
    frontTime = frontTimePix/BigPicRotStruct.fps;       % time vector for front
    [~,fronLocOnTime] = min(abs(bsxfun(@minus,BigPicRotStruct.t(:),frontTime(:)')));
    %-- front's bottom:
    bottomDistPix = round(bottomDist./BigPicRotStruct.res);
    dxPix = round(dx./BigPicRotStruct.res);
    regionsForBottomVal_x =  bsxfun(@plus,repmat(((bottomDistPix-dxPix):(bottomDistPix+dxPix))+1,length(BigPicRotStruct.x(:)),1),(0:(length(BigPicRotStruct.x)-1))');
    regionsForBottomVal_x(regionsForBottomVal_x<1) = 1;
    regionsForBottomVal_x(regionsForBottomVal_x>length(BigPicRotStruct.x)) = length(BigPicRotStruct.x);
    regionsForBottomVal_t = round(repmat(fronLocOnTime(:),1,size(regionsForBottomVal_x,2)));
    linearInd = sub2ind(size(BigPicRotStruct.DataMatNorm), regionsForBottomVal_t, regionsForBottomVal_x);
    
    %--smooth
    DataMatSmth = conv2(BigPicRotStruct.DataMatNorm,ones(1,smt)/smt,'same');
    
    bottomVal = mean(DataMatSmth(linearInd),2);
    
    %--- create normalized and streched matrix
    DataMatShifted = bsxfun(@minus,DataMatSmth,bottomVal(:)');
    DataMatStreched = bsxfun(@rdivide,DataMatShifted,1-bottomVal(:)');
    Xc_vec = nan(size(BigPicRotStruct.x));
    for i=XX
        Tidx = fronLocOnTime(i);
        [Xi,~] = intersections(XX,DataMatStreched(Tidx,:),[-inf,inf],Afall*[1 1]);
        %-- keep only crossing location between front and 'bottomDist'
        Xkeep = Xi(Xi<=i & Xi>=(i+bottomDistPix));
        if ~isempty(Xkeep)
            %-- keep closest to x_tip:
            Xc_vec(i) = (i-max(Xkeep))*BigPicRotStruct.res;
        else
            continue
        end
    end
    
    
    chsvFromContact.Xc = Xc_vec;
    chsvFromContact.DA = 1-bottomVal;
    chsvFromContact.Afall = Afall;
    chsvFromContact.bottomDist = bottomDist;
    chsvFromContact.dx = dx;
    
%     %% by single frame
%     %-- find frame
%     [X0,T0] = intersections(BigPicRotStruct.x,frontTime,PhotoLocation*[1 1],[-1,1]);    % get time when from reches 'PhotoLocation'
%     [~,I] = min(abs(BigPicRotStruct.t-T0));
%     X_minXtip = BigPicRotStruct.x-X0;
%     snapsht = BigPicRotStruct.DataMatNorm(I,:);   %snapshot relevant to requested location 'PhotoLocation'
%     [~,BottmPointPix] = min(abs(X_minXtip-bottomDist));
%     deltaM = 1e-3; % define the width of window on which to measure the bottom value
%     deltaP = round(deltaM/BigPicRotStruct.res);
%     BottmPointVal = mean(snapsht(BottmPointPix-deltaP:BottmPointPix+deltaP));
%     snapshtStrchd = (snapsht-BottmPointVal)./(1-BottmPointVal);
%     snapshtSmth = smooth(snapshtStrchd,smoothMask);
%     [XcPosble,~] = intersections(X_minXtip,snapshtSmth,[-100,100],Afall*[1 1]);
%     Xc = max(XcPosble(XcPosble<0));
%     if isempty(Xc)
%         Xc = nan;
%     end
%     
%     
%     [~,iLoc] = min(abs(smooth(BigPicRotStruct.frontVelLoc_interpM,5)-PhotoLocation));
%     Cf = BigPicRotStruct.frontVel_interpMperS(iLoc);
%     
%     Tc = -Cf.*Xc;
%     
%     [~,zoerPointPix] = min(abs(X_minXtip+bottomDist));
%     A0 = mean(BigPicRotStruct.DataMat(I,zoerPointPix-deltaP:zoerPointPix+deltaP));
%     DA = A0 - mean(BigPicRotStruct.DataMat(I,BottmPointPix-deltaP:BottmPointPix+deltaP));
%     
%     chsvFromContact.Xc = Xc;
%     chsvFromContact.Tc = Tc;
%     chsvFromContact.DA = DA;
%     chsvFromContact.A0 = A0;
%     chsvFromContact.Cf = Cf;
%     chsvFromContact.Yfall = Afall;
%     chsvFromContact.bottomLoc = bottomDist;
%     chsvFromContact.rotPhotoRowNorm = snapshtSmth;
    
end