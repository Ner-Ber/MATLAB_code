function phedi_findCohesiveZoneKink(PhediData,varargin)


%% set defaults
[KinkSection] = setDefaults4function(varargin,[-4 -1]*1e-4);

%% get relevant section of kink
t_mins_t_tips = bsxfun(@minus,PhediData.timeVec,PhediData.t_tips(:)');
N = size(t_mins_t_tips,2);      % number of phedis
KinkTime = t_mins_t_tips>=KinkSection(1) &t_mins_t_tips<=KinkSection(2);
t_mins_t_tipsLinear = t_mins_t_tips(KinkTime);
t_mins_t_tipsSection = reshape(t_mins_t_tipsLinear,[],N);
PhediLocZeros = bsxfun(@minus,PhediData.PhediLocation,mean(PhediData.PhediLocation(1:100,:),1,'omitnan'));
PhediLocZerosLinear = PhediLocZeros(KinkTime);
PhediLocSection = reshape(PhediLocZerosLinear,[],N);

%% shift the signal to origin of axis
timeShifts = t_mins_t_tipsSection(1,:);
t_mins_t_tipsShifted = bsxfun(@minus,t_mins_t_tipsSection,timeShifts);
t_strech = t_mins_t_tipsShifted(end,:);
t_mins_t_tipsStreched = bsxfun(@rdivide,t_mins_t_tipsShifted,t_strech);

ampShifts = PhediLocSection(1,:);
PhediLocationShifted = bsxfun(@minus,PhediLocSection,ampShifts);
ampStrech = PhediLocationShifted(end,:);
PhediLocationStreched = bsxfun(@rdivide,PhediLocationShifted,ampStrech);

L = size(PhediLocationStreched,1);

%% iterate on phedis to find fits
%-- prepare data for optimization
% MaxAmp = max(PhediLocationShifted(:));
% MinAmp = min(PhediLocationShifted(:));
% MaxAmpMean = mean(max(PhediLocationShifted));
% MinAmpMean = mean(min(PhediLocationShifted));
% maxTime = max(t_mins_t_tipsShifted(:));
% p2pSlope = (MaxAmpMean-MinAmpMean)./maxTime;
% P0 = [0,0, p2pSlope, maxTime/2];
% epsilon = maxTime/100;
% ub = [MaxAmp, p2pSlope*10, p2pSlope*10, maxTime-epsilon];
% lb = [MinAmp, -p2pSlope*10, -p2pSlope*10, epsilon];

kinkLocVec = nan(N,1);
P_cell = cell(N,1);
for i = 1:N
    %--- set initial and boudaries
%     MaxAmp = max(PhediLocationStreched(:,i));
%     MinAmp = min(PhediLocationStreched(:,i));
%     maxTime = max(t_mins_t_tipsStreched(:,i));
%     maxSlope_r = (MaxAmp-MinAmp)/(0.1*maxTime);
%     minSlope_r = 0;
%     initSlope_r = (MaxAmp-PhediLocationStreched(round(L/2),i))/(maxTime/2);
%     maxSlope_l = maxSlope_r;
%     minSlope_l = -maxSlope_r;
%     initSlope_l = mean([MaxAmp,MinAmp])./(maxTime/2);
%     initialKinkLoc = maxTime/2;
%     initialShift = PhediLocationStreched(round(L/2),i);
%     epsilon = maxTime/100;
    
    %     P0 = [initialShift, initSlope_l, initSlope_r, initialKinkLoc];
    %     ub = [MaxAmp, maxSlope_l, maxSlope_r, maxTime-epsilon];
    %     lb = [MinAmp, minSlope_l, minSlope_r, epsilon];
    
%     P0 = [initSlope_l, initSlope_r, initialKinkLoc];
%     ub = [maxSlope_l, maxSlope_r, maxTime-epsilon];
%     lb = [minSlope_l, minSlope_r, epsilon];

    P0 = [0.5, 0];
    ub = [0.95, 1];
    lb = [0.05, -1];

    plusfun = @(x) max(x,0);
%     model = @(P,x) P(1)*plusfun(P(3)-x) + P(2)*plusfun(x-P(3)) + P(1)*abs(P(3));
    model = @(P,x) (P(2)/P(1))*plusfun(P(1)-x) + ((1-P(2))./(1-P(1)))*plusfun(x-P(1)) + P(2);
    
    P = lsqcurvefit(model,P0,t_mins_t_tipsStreched(:,i),PhediLocationStreched(:,i),...
        lb,ub);
    
    %     P = fit_piecewiseLinear(t_mins_t_tipsShifted(:,i),PhediLocationShifted(:,i),...
    %         P0,lb,ub);
    kinkLocVec(i) = P(1);
    P_cell{i} = P;
end

end