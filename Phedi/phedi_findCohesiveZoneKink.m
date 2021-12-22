function [X_kink,Y_kink] = phedi_findCohesiveZoneKink(PhediData,varargin)
% [X_kink,Y_kink] = phedi_findCohesiveZoneKink(PhediData,KinkSection)
%
% KinkSection - the time region in which to search for the kink ([-4
% -0.2]*1e-4 is the defaults).


%% set defaults
[KinkSection] = setDefaults4function(varargin,[-4 -0.2]*1e-4);

%% get relevant section of kink
t_mins_t_tips = bsxfun(@minus,PhediData.timeVec(:),PhediData.t_tips(:)');
N = size(t_mins_t_tips,2);      % number of phedis
KinkTime = t_mins_t_tips>=KinkSection(1) &t_mins_t_tips<=KinkSection(2);

[~,col] = find(KinkTime);
k = find(KinkTime);
[C] = unique(col);
Ccell = num2cell(C);
K = cellfun(@(A) k(A==col),Ccell,'UniformOutput',0);
t_mins_t_tipsSection = cellfun(@(A) t_mins_t_tips(A),K,'UniformOutput',0);

PhediLocZeros = bsxfun(@minus,PhediData.PhediLocation,mean(PhediData.PhediLocation(1:100,:),1,'omitnan'));
% PhediLocZerosLinear = PhediLocZeros(KinkTime);

PhediLocSection = cellfun(@(A) PhediLocZeros(A),K,'UniformOutput',0);

% PhediLocSection = reshape(PhediLocZerosLinear,[],N);

%% shift the signal to origin of axis and strech to normalize
% timeShifts = mean(t_mins_t_tipsSection(1:20,:));
timeShifts = cellfun(@(A) A(1,:),t_mins_t_tipsSection);
% t_mins_t_tipsShifted = bsxfun(@minus,t_mins_t_tipsSection,timeShifts);
t_mins_t_tipsShifted = cellfun(@(A,B) A-B,t_mins_t_tipsSection,num2cell(timeShifts),'UniformOutput',0);
% t_strech = mean(t_mins_t_tipsShifted(end-20:end,:));
t_strech = cellfun(@(A) A(end),t_mins_t_tipsShifted);

t_mins_t_tipsStreched = cellfun(@(A) A./A(end),t_mins_t_tipsShifted,'UniformOutput',0);

ampShifts  = cellfun(@(A) A(1),PhediLocSection);
PhediLocationShifted = cellfun(@(A) A-A(1),PhediLocSection,'UniformOutput',0);
ampStrech = cellfun(@(A) A(end),PhediLocationShifted);
PhediLocationStreched = cellfun(@(A) A/A(end),PhediLocationShifted,'UniformOutput',0);


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
    P = lsqcurvefit(model,P0,t_mins_t_tipsStreched{i},PhediLocationStreched{i},...
        lb,ub,opts);
    
    X_kink_norm(i) = P(1);
    Y_kink_norm(i) = P(2);
end

%% restrech and renormalize kink location
X_kink = X_kink_norm(:).*t_strech(:)+timeShifts(:);
Y_kink = Y_kink_norm(:).*ampStrech(:)+ampShifts(:);


end