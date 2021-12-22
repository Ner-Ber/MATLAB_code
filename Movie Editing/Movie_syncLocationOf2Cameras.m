function [bigPicLocationPix, RMS_vec] = Movie_syncLocationOf2Cameras(t_big,ROT_big,t_small,ROT_small,varargin)
    % [bigPicLocationPix, RMS_vec] = Movie_syncLocationOf2Cameras(t_big,ROT_big,t_small,ROT_small,varargin)
    %varargin: RotBigAvgWindow,TpostFall,Nregions] = [11,5e-5,1]
    %
    % will find the picel in the BigPic of where the center of the SmallPic
    % (phedis photograph) is located. This is donce by fitting the falling
    % of contact area.
    
    %% set defaults
    [RotBigAvgWindow,TpostFall,Nregions] = setDefaults4function(varargin,11,5e-5,1);
    
    %% normalize to prepare for fit
    %--- Norm BigPic:
    ROT_bigNorm = bsxfun(@rdivide,ROT_big,mean(ROT_big(1:100,:)));
    ROT_bigNormAvg = movmean(ROT_bigNorm,RotBigAvgWindow,2);
    %--- Norm SmallPic:
    Lregion = floor(size(ROT_small,2)/Nregions);
    D = Lregion*(1:Nregions)';
    DD = [D-Lregion+1,D];
    DD(end,end) = size(ROT_small,2);
    ROT_smallAvgVec = nan(size(ROT_small,1),Nregions);
    ROT_smallAvgVecNrom = nan(size(ROT_small,1),Nregions);
    ROT_smallStrch = nan(size(ROT_small,1),Nregions);
    Tinter = zeros(Nregions,1);
    bigPicLocationPix_vec = nan(Nregions,1);
    
    %--- find falling edge region of the front:
    for i=1:Nregions
         R_temp = mean(ROT_small(:,[DD(i,1):DD(i,2)]),2);
        R_temp2 = R_temp./mean(R_temp(1:100));
        Iinter = find(R_temp2<0.95,1,'first');
        if isempty(Iinter)
            warning('No falling edge found in small ROT. Might be precoursor');
            RMS_vec = nan;
            break
        end
        Tinter(i) = t_small(Iinter);
        
        %--- normalize before fall:
        reg2normLogic = t_small<(Tinter(i)-TpostFall) & t_small>(Tinter(i)-2*TpostFall);
        ROT_smallAvgVecNrom(:,i) = R_temp./mean(R_temp(reg2normLogic));
        
        
        %--- stretch from 0 to 1:
        [~,Ibig] = min(abs(t_big-(Tinter(i)+TpostFall)));
        ROT_bigStrch = bsxfun(@rdivide,bsxfun(@minus, ROT_bigNormAvg,ROT_bigNormAvg(Ibig,:)),(1-ROT_bigNormAvg(Ibig,:)));
        [~,Ismall] = min(abs(t_small-(Tinter(i)+TpostFall)));
        ROT_smallStrch(:,i) = (ROT_smallAvgVecNrom(:,i)-ROT_smallAvgVecNrom(Ismall,i))./(1-ROT_smallAvgVecNrom(Ismall,i));
        
        
        %% find column with best fit
        %--- create common time axis:
        timeRegion = [Tinter(i)-TpostFall,Tinter(i)+TpostFall];
        timeRegion_LogicSmall = t_small>=min(timeRegion) & t_small<=max(timeRegion);
        timeRegion_LogicBig = t_big>=min(timeRegion) & t_big<=max(timeRegion);
        
        unifiedStart = max(min(t_small(timeRegion_LogicSmall)),min(t_big(timeRegion_LogicBig)));
        unifiedEnd = min(max(t_small(timeRegion_LogicSmall)),max(t_big(timeRegion_LogicBig)));
        unifiedTime = linspace(unifiedStart,unifiedEnd,2*max(nnz(timeRegion_LogicSmall),nnz(timeRegion_LogicBig)));
        A_small = interp1(t_small,ROT_smallStrch(:,i),unifiedTime)';
        [Xo,To] = meshgrid(1:size(ROT_bigStrch,2),t_big(timeRegion_LogicBig));
        [Xq,Tq] = meshgrid(1:size(ROT_bigStrch,2),unifiedTime);
        A_big = interp2(Xo,To,ROT_bigStrch(timeRegion_LogicBig,:),Xq,Tq);

        %-- find minimal rms:
        subtracted = bsxfun(@minus,A_big,A_small);
        [RMSi,minRMSpix] = min(rms(subtracted,1));
        RMS_vec(i) = RMSi;
        bigPicLocationPix_vec(i) = minRMSpix;

    end
    
    bigPicLocationPix = mean(bigPicLocationPix_vec);
end