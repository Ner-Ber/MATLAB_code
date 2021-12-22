function [X_kink,Y_kink] = phedi_findCohesiveZoneKink(PhediData,varargin)


%% set defaults
[KinkSection] = setDefaults4function(varargin,[-4 -0.2]*1e-4);

%% get relevant section of kink
t_mins_t_tips = bsxfun(@minus,PhediData.timeVec,PhediData.t_tips(:)');
N = size(t_mins_t_tips,2);      % number of phedis
KinkTime = t_mins_t_tips>=KinkSection(1) &t_mins_t_tips<=KinkSection(2);
t_mins_t_tipsLinear = t_mins_t_tips(KinkTime);
t_mins_t_tipsSection = reshape(t_mins_t_tipsLinear,[],N);
PhediLocZeros = bsxfun(@minus,PhediData.PhediLocation,mean(PhediData.PhediLocation(1:100,:),1,'omitnan'));
PhediLocZerosLinear = PhediLocZeros(KinkTime);
PhediLocSection = reshape(PhediLocZerosLinear,[],N);

%% shift the signal to origin of axis and strech to normalize
% timeShifts = mean(t_mins_t_tipsSection(1:20,:));
timeShifts = t_mins_t_tipsSection(1,:);
t_mins_t_tipsShifted = bsxfun(@minus,t_mins_t_tipsSection,timeShifts);
% t_strech = mean(t_mins_t_tipsShifted(end-20:end,:));
t_strech = t_mins_t_tipsShifted(end,:);
t_mins_t_tipsStreched = bsxfun(@rdivide,t_mins_t_tipsShifted,t_strech);

% ampShifts = mean(PhediLocSection(1:20,:));
ampShifts = PhediLocSection(1,:);
PhediLocationShifted = bsxfun(@minus,PhediLocSection,ampShifts);
% ampStrech = mean(PhediLocationShifted(end-20:end,:));
ampStrech = PhediLocationShifted(end,:);
PhediLocationStreched = bsxfun(@rdivide,PhediLocationShifted,ampStrech);


%% iterate on phedis to find fits
%--- set initial and boudaries
P0 = [0.5, 0];
ub = [0.95, 1];
lb = [0.05, -1];

X_kink_norm = nan(N,1);
Y_kink_norm = nan(N,1);
for i = 1:N
    model = @(P,x) heaviside(P(1)-x).*(P(2)/P(1)).*x + heaviside(x-P(1)).*(((1-P(2)).*x+P(2)-P(1))./(1-P(1)));
    
    opts = optimset('Display','off');
    P = lsqcurvefit(model,P0,t_mins_t_tipsStreched(:,i),PhediLocationStreched(:,i),...
        lb,ub,opts);
    
    X_kink_norm(i) = P(1);
    Y_kink_norm(i) = P(2);
end

%% restrech and renormalize kink location
X_kink = X_kink_norm(:).*t_strech(:)+timeShifts(:);
Y_kink = Y_kink_norm(:).*ampStrech(:)+ampShifts(:);


end