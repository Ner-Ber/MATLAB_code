function [P] = fit_piecewiseLinear(xData,yData,varargin)
% [P] = FIT_PIECEWISELINEAR (xData,yData,P0, lb, ub)
%FIT_PIECEWISELINEAR will fit a model of two-linear-parts to the given x,y
%data.
%
% the model is given by
% model = @(P,x) P(1) + P(2)*plusfun(P(4)-x) + P(3)*plusfun(x-P(4))
% where plusfun=max(x,0);
% xData and yData are mandatory inputs. all other are optional. 
% P0 - initial values for model fitting
% lb - lower bounds
% ub - upper boubds

P0Def = [0 0 0.003 1.5e-4];
lbDef = [-0.1 -0.5 0 0];
ubDef = [0.1 0.5 0.5 0.1];

[P0, lb, ub] = setDefaults4function(varargin,P0Def, lbDef, ubDef);

plusfun = @(x) max(x,0);
model = @(P,x) P(1) + P(2)*plusfun(P(4)-x) + P(3)*plusfun(x-P(4));
P = lsqcurvefit(model,P0,xData,yData,lb,ub);


end