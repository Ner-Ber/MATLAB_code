function [XportCoeffs, fresnelComplex] = phedi_fitFresnelFunction(x,y,varargin)
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

%-- peak location of unshifted unsreched function:
[YPeak_fun, XPeak_fun] = deal(1.370442919651484,-0.971180180180180);
[YLowPoint_fun, XLowPoint_fun] = deal(0.007967611937845, 2);

%-- get initial conditions for fit:
[Ypeak_data, loc_p] = max(y);
XPeak_data = x(loc_p);
[YLowPoint_data, loc_l] = min(y);
XLowPoint_data = x(loc_l);

InitialyStrech = Ypeak_data/YPeak_fun;
InitialxStrech = (XLowPoint_fun-XPeak_fun)./(XLowPoint_data-XPeak_data);
initialxDisplcmnt = XPeak_data - XPeak_fun/InitialxStrech;
initialyDisplcmnt = 0;

x0 = [InitialyStrech; InitialxStrech; initialxDisplcmnt; initialyDisplcmnt];

% fresnelComplex = @(b,t)...
%     b(1)*(abs((1+1i)./(2^1.5)-(2^-0.5).*(fresnelc(sqrt(pi/2)*b(2)*(t-b(3)))+1i*fresnels(sqrt(pi/2)*b(2)*(t-b(3)))))).^2 +b(4);
fresnelComplex = @(xo,xdata) fitFresnelIntensity(xo,xdata);
%% set given parameters
lb = [vaerticalStrech, horizontalStrech, horizontalDisplacement,verticalDisplacement];
ub = [vaerticalStrech, horizontalStrech, horizontalDisplacement,verticalDisplacement];
lb(isnan(lb)) = -inf;
ub(isnan(ub)) = +inf;

%% find fit
options = optimoptions('lsqcurvefit','Display','none');
XportCoeffs = lsqcurvefit(fresnelComplex,x0,x,y,lb,ub,options);
% mdl = fitnlm(x,y,fresnelComplex,x0);
% XportCoeffs = mdl.Coefficients{:,1};

end