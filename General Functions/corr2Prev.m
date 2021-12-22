function [lags, corr_val, N_vector] = corr2Prev(vectorMatrix, n)
% [lags, corr_val, N_vector] = corr2Prev(vectorMatrix)
% 
% corr2Prev will compare each column in vectorMatrix to the n previous
% column to find the maximal correlation among them n. In case there are
% less than n previous columns the current column will be compared to all
% the previous ones. 
% all lages are given in comparison to the first column
%
% vectorMatrix - a MxN matric, each column is a time series composed of M
% samples. 
% n - (optional) number of previous columns to compare to. Default is n=2.
% n can be stated as n='all' in which case n=N. 
% 

%% set defaults
[~, N] = size(vectorMatrix);
vectorMatrix(isnan(vectorMatrix)) = 0;
if nargin<2
    n=2;
elseif ischar(n)
    n=N;
end

%% uterate upon columns
[lags, corr_val, N_vector] = deal(0);

for i = 2:N
    ni = min(n,i-1);
    %--- iterate upon all previous columns
    RK = -inf(ni,1);
    LAGK = -nan(ni,1);
    for k = 1:ni
        [rk,lag_k] = xcorr(vectorMatrix(:,i),vectorMatrix(:,i-k));
        [rk_max, rk_m_idx] = max(rk);
        RK(k) = rk_max;
        LAGK(k) = lag_k(rk_m_idx);
    end
    
    [ri,IDXi] = max(RK);
    lags(i) = LAGK(IDXi);
    corr_val(i) = ri;
    N_vector(i) = i-IDXi;
    
end


end