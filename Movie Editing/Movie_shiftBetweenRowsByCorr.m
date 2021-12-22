function [S,varargout] = Movie_shiftBetweenRowsByCorr(RowOverTime,varargin)
%Movie_shiftBetweenRowsByCorr is meant to find the (approx.) frame of
%shifting by correlation between rows.
% OPTIONAL: row to consider.

[M,~] = size(RowOverTime);
[RelevantRow] = setDefaults4function(varargin,1:M);
ROT_red = RowOverTime(RelevantRow,:);

[N,~] = size(ROT_red);
L = zeros(N,0);
% R = L;
for i = 2:N
    [r,lags] = xcorr(ROT_red(1,:),ROT_red(i,:));
    [~,I] = max(r);
%     R(i) = Ir;
    L(i) = lags(I);
end

Si = find(~L,1,'last');
varargout{1} = Si;
S = RelevantRow(Si);