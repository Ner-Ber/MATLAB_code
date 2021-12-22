function [parameters] = analyzeFindSqrtDecayParmtrs(x_ax,y_ax)
%y_ax - contains from peak forward

LogLikeFunc = @(x,a,b,c) sum(log((a./sqrt(x+b))+c));

