function [power,sqrtPref] = Crack_UxFromInterfaceSolution(itnerfaceSolutionStruct)

Ux_fromvx = -(1/itnerfaceSolutionStruct.v)*cumtrapz(itnerfaceSolutionStruct.x,itnerfaceSolutionStruct.vx);

%--- log log fit
Ux_fromvxShifted = Ux_fromvx-min(Ux_fromvx);
positiveLogical = itnerfaceSolutionStruct.x<0;
XX = -itnerfaceSolutionStruct.x(positiveLogical);
UxUx = Ux_fromvxShifted(positiveLogical);
%--- take largest half of X valuews
[~,I] = sort(XX);
XXlarge = XX(I(round(length(I)/2):end));
UxUxLarge = UxUx(I(round(length(I)/2):end));

Plarge = polyfit(log10(XXlarge),log10(UxUxLarge),1);
power = Plarge(1);
sqrtPref = 10^Plarge(2);
end