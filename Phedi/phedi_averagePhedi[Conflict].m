function AvgPhediStruct = phedi_averagePhedi(DataStruct, varargin)
    % AvgPhediStruct = phedi_averagePhedi(DataStruct, phediIdxs, reduceSin, t_tips, plotSinFit)
    %
    %INPUTS:
    % DataStruct - the main output of phedi_createStructureWAdd.
    % phediIdxs - indexes of phedis to calculate their common value (default=all).
    %
    %OPTIONAL:
    % reduceSin, t_tips, plotSinFits
    
    %% get phedi data
    PhediVel = DataStruct.PhediData.PhediVelocity;
    N = size(PhediVel,2);
    PhediLoc = DataStruct.PhediData.PhediLocation; 
    
    %% set defaults
    [phediIdxs, reduceSin, t_tips, plotSinFit] = setDefaults4function(varargin,true(1,N),1,1,0);
    
    %% create fixed Time matrix for
    if t_tips
        T1 = DataStruct.PhediData.t_mins_t_tip_4vel;
        T2 = DataStruct.PhediData.t_mins_t_tip;
        
        X1 = DataStruct.PhediData.x_mins_x_tip_4vel;
        X2 = DataStruct.PhediData.x_mins_x_tip;
    else
        T1 = repmat(DataStruct.PhediData.timeVec_4vel,1,N);
        T2 = repmat(DataStruct.PhediData.timeVec,1,N);
        
        if length(T1)~=length(T2)
            X1 = repmat(movmean(DataStruct.PhediData.PhediLocation,2,'Endpoints','discard'),1,N);
        else
            X1 = repmat(DataStruct.PhediData.PhediLocation,1,N);
        end
        X2 = repmat(DataStruct.PhediData.PhediLocation,1,N);
    end
    %-- eliminate nans
    X1(abs(X1)==inf) = nan;
    X2(abs(X2)==inf) = nan;
    T1(abs(T1)==inf) = nan;
    T2(abs(T2)==inf) = nan;
%     M = size(X1,1);
%     [row,~] = find(abs(X1)==inf);
%     Min = max(row(row<=M/2));   if isempty(Min); Min=0; end;
%     Max = min(row(row>=M/2));   if isempty(Max); Max=M; end;
%     X1 = X1((Min+1):(Max-1),:);
%     M = size(X2,1);
%     [row,~] = find(abs(X2)==inf);
%     Min = max(row(row<=M/2));   if isempty(Min); Min=0; end;
%     Max = min(row(row>=M/2));   if isempty(Max); Max=M; end;
%     X2 = X2((Min+1):(Max-1),:);
%     M = size(T1,1);
%     [row,~] = find(abs(T1)==inf);
%     Min = max(row(row<=M/2));   if isempty(Min); Min=0; end;
%     Max = min(row(row>=M/2));   if isempty(Max); Max=M; end;
%     T1 = T1((Min+1):(Max-1),:);
%     M = size(T2,1);
%     [row,~] = find(abs(T2)==inf);
%     Min = max(row(row<=M/2));   if isempty(Min); Min=0; end;
%     Max = min(row(row>=M/2));   if isempty(Max); Max=M; end;
%     T2 = T2((Min+1):(Max-1),:);
    
    T = sort([T1;T2]);
    X = sort([X1;X2]);
    
     %--- reduce sin if needed
    if reduceSin
        myfitStruct = phedi_fitSineToPhediTraj(DataStruct.PhediData,[],2,plotSinFit);
        PhediLocFixed = nan(size(PhediLoc));
        for i=1:N
            PhediLocFixed(:,i) = PhediLoc(:,i)-myfitStruct.model(myfitStruct.b_param{i},myfitStruct.Shifts(i),T2(:,i));
        end
    else
        PhediLocFixed = PhediLoc;
    end
    
    %% create common data
    
    %--- interpolate for common time steps
    t_unique = unique(T);
    T_vec = min(t_unique):abs(min(min(diff(T2))))/4:max(t_unique);
    T_phedis = repmat(T_vec(:),1,N);
    %--- interpolate for common spatial steps
    x_unique = unique(X);
    %-- exclude nans:
    DD = diff(X2);
%     X_vec = min(x_unique(abs(x_unique)~=inf)):(abs(min(min(DD(abs(DD)~=inf))))/4):max(x_unique(abs(x_unique)~=inf));
    minBound = max(min(x_unique(abs(x_unique)~=inf)),-10e-2);
    maxBound = min(max(x_unique(abs(x_unique)~=inf)),10e-2);
%     IntrpStep = max((min(min(abs(DD(abs(DD)~=inf))))/4),5e-6);
    IntrpStep = max((min(min(abs(DD(abs(DD)~=inf))))/4),5e-4);
    X_vec = minBound:IntrpStep:maxBound;
    X_phedis = repmat(X_vec(:),1,N);

    
    
    PhediLoc_iterped_t = nan(size(T_phedis));
    PhediVel_iterped_t = nan(size(T_phedis));
    PhediLoc_iterped_x = nan(size(X_phedis));
    PhediVel_iterped_x = nan(size(X_phedis));
    %--- interp location and velocity for time and spacr
    for i=1:N
        NanIdI = isnan(T2(:,i)) | isinf(T2(:,i));
        XI = T2(~NanIdI,i); VI = PhediLocFixed(~NanIdI,i);
        NanIdQ = isnan(T_phedis(:,i)) | isinf(T_phedis(:,i));
        XQ = T_phedis(:,i); XQ(NanIdQ) = 0;
        VQ = interp1(XI, VI,  XQ);
        VQ(NanIdQ) = nan; PhediLoc_iterped_t(:,i) = VQ;
        
        NanIdI = isnan(T1(:,i)) | isinf(T1(:,i));
        XI = T1(~NanIdI,i); VI = PhediVel(~NanIdI,i);
        NanIdQ = isnan(T_phedis(:,i)) | isinf(T_phedis(:,i));
        XQ = T_phedis(:,i); XQ(NanIdQ) = 0;
        VQ = interp1(XI, VI,  XQ);
        VQ(NanIdQ) = nan; PhediVel_iterped_t(:,i) = VQ;
        
        NanIdI = isnan(X2(:,i)) | isinf(X2(:,i));
        XI = X2(~NanIdI,i); VI = PhediLocFixed(~NanIdI,i);
        NanIdQ = isnan(X_phedis(:,i)) | isinf(X_phedis(:,i));
        XQ = X_phedis(:,i); XQ(NanIdQ) = 0;
        VQ = interp1(XI, VI,  XQ);
        VQ(NanIdQ) = nan; PhediLoc_iterped_x(:,i) = VQ;
        
        NanIdI = isnan(X1(:,i)) | isinf(X1(:,i));
        XI = X1(~NanIdI,i); VI = PhediVel(~NanIdI,i);
        NanIdQ = isnan(X_phedis(:,i)) | isinf(X_phedis(:,i));
        XQ = X_phedis(:,i); XQ(NanIdQ) = 0;
        VQ = interp1(XI, VI,  XQ);
        VQ(NanIdQ) = nan; PhediVel_iterped_x(:,i) = VQ;
        
%         PhediLoc_iterped_t(:,i) = interp1(T2(:,i), PhediLocFixed(:,i),  T_phedis(:,i));
%         PhediVel_iterped_t(:,i) = interp1(T1(:,i), PhediVel(:,i),       T_phedis(:,i));
%         PhediLoc_iterped_x(:,i) = interp1(X2(:,i), PhediLocFixed(:,i),  X_phedis(:,i));
%         PhediVel_iterped_x(:,i) = interp1(X1(:,i), PhediVel(:,i),       X_phedis(:,i));
    end
    
    %--- create measurements for time
    T_mean = mean(T_phedis(:,phediIdxs),2);
    PhediLoc_mean_t = mean(PhediLoc_iterped_t(:,phediIdxs),2);
    PhediLoc_var_t = var(PhediLoc_iterped_t(:,phediIdxs),0,2);
    PhediLoc_max_t = max(PhediLoc_iterped_t(:,phediIdxs),[],2);
    PhediLoc_min_t = min(PhediLoc_iterped_t(:,phediIdxs),[],2);
    PhediVel_mean_t = mean(PhediVel_iterped_t(:,phediIdxs),2);
    PhediVel_var_t = var(PhediVel_iterped_t(:,phediIdxs),0,2);
    PhediVel_max_t = max(PhediVel_iterped_t(:,phediIdxs),[],2);
    PhediVel_min_t = min(PhediVel_iterped_t(:,phediIdxs),[],2);
    
    %--- create measurements for space
    X_mean = mean(X_phedis(:,phediIdxs),2);
    PhediLoc_mean_x = mean(PhediLoc_iterped_x(:,phediIdxs),2);
    PhediLoc_var_x = var(PhediLoc_iterped_x(:,phediIdxs),0,2);
    PhediLoc_max_x = max(PhediLoc_iterped_x(:,phediIdxs),[],2);
    PhediLoc_min_x = min(PhediLoc_iterped_x(:,phediIdxs),[],2);
    PhediVel_mean_x = mean(PhediVel_iterped_x(:,phediIdxs),2);
    PhediVel_var_x = var(PhediVel_iterped_x(:,phediIdxs),0,2);
    PhediVel_max_x = max(PhediVel_iterped_x(:,phediIdxs),[],2);
    PhediVel_min_x = min(PhediVel_iterped_x(:,phediIdxs),[],2);
    
    %% save to structure
    AvgPhediStruct = struct(...
            'T_mean',T_mean,...
            'Loc_mean_t', PhediLoc_mean_t,...
            'Loc_var_t',PhediLoc_var_t,...
            'Loc_max_t', PhediLoc_max_t,...
            'Loc_min_t',PhediLoc_min_t,...
            'Vel_mean_t',PhediVel_mean_t,...
            'Vel_var_t',PhediVel_var_t,...
            'Vel_max_t',PhediVel_max_t,...
            'Vel_min_t',PhediVel_min_t,...
            ...
            'X_mean',X_mean,...
            'Loc_mean_x', PhediLoc_mean_x,...
            'Loc_var_x',PhediLoc_var_x,...
            'Loc_max_x', PhediLoc_max_x,...
            'Loc_min_x',PhediLoc_min_x,...
            'Vel_mean_x',PhediVel_mean_x,...
            'Vel_var_x',PhediVel_var_x,...
            'Vel_max_x',PhediVel_max_x,...
            'Vel_min_x',PhediVel_min_x...
            );
end