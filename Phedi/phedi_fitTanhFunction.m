function [XportCoeffs, tanhFit] = phedi_fitTanhFunction(x,y,varargin)
%[XportCoeffs, fresnelComplex] = phedi_fitFresnelFunction(x,y,vaerticalStrech, horizontalStrech, horizontalDisplacement,verticalDisplacement)
% phedi_fitFresnelFunction will fit a complex Frensel function to a given
% data.
%
%   inputs:
%   - x,y: coordinates of the data to fit to.
%
%   outputs:
%   -XportCoeffs: coefficients of the fit by the order:
%       vaertical strech, horizontal strech, horizontal displacement,
%       vertical displacement
%   -fresnelComplex: function handle to recieve XportCoeffs, in the format:
%       fresnelComplex = @(coeffs,x)...

%% set given parameters int vector

[vaerticalStrech, horizontalStrech, horizontalDisplacement,verticalDisplacement] = setDefaults4function(varargin,nan,nan,nan,nan);

%%


%-- get initial conditions for fit:
InitialyStrech = max(y)-min(y);
InitialxStrech = 1;
initialxDisplcmnt = mean(x);
initialyDisplcmnt = mean(y);


x0 = [InitialyStrech; InitialxStrech; initialxDisplcmnt; initialyDisplcmnt];

tanhFit = @(b,t)...
b(1)*tanh(b(2)*(t-b(3)))+b(4);
%% set given parameters
lb = [vaerticalStrech, horizontalStrech, horizontalDisplacement,verticalDisplacement];
ub = [vaerticalStrech, horizontalStrech, horizontalDisplacement,verticalDisplacement];
lb(isnan(lb)) = -inf;
ub(isnan(ub)) = +inf;

%% find fit
options = optimoptions('lsqcurvefit','Display','none');
XportCoeffs = lsqcurvefit(tanhFit,x0,x,y,lb,ub,options);

end