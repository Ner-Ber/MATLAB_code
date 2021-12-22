function [selectedShift, newX] = phedi_fitToPowerAndAmp(x,y,power,amplitude,varargin)
%[fitPowerStruct] = phedi_fitToPowerAndAmp(x,y,power,amplitude, x_min,x_max, Shifts)
%
% phedi_fitToPowerAndAmp is used to find the shift s in the expression
% y=A*(x-s)^p. Here A and p are the amplitude and power which are required
% for this function.
%
% DAEFAULTS
% x_min = 1e-5;
% x_max = 9e-5;
% Shifts = linspace(-2*1e-3,2*1e-3,150);


%% set defaults
warning('off');
%% set data
%-- parameters:
x_min_def = 1e-5;
x_max_def = 9e-5;
Shifts_def = linspace(-2*1e-3,2*1e-3,1e4);
[x_min, x_max, Shifts] = setDefaults4function(varargin, x_min_def,x_max_def, Shifts_def);

%%
xLogical = x<=x_max & x>=x_min;
XX = x(xLogical);
YY = y(xLogical);
RMS_vec = nan(1,length(Shifts));
for j = 1:length(Shifts)
    s = Shifts(j);
    shiftedXX = XX-s;
    positiveLogic = shiftedXX>0;
    model = amplitude*((shiftedXX(positiveLogic)).^power);
    measure = YY(positiveLogic);
    RMS_vec(j) = rms(model - measure);
end

%% find minimal RMS:
[~,I] = min(RMS_vec);
if I==1 || I==length(Shifts)
    selectedShift = Shifts(I);
else
    P = polyfit(Shifts(I-1:I+1),RMS_vec(I-1:I+1),2);
    selectedShift = -P(2)/(2*P(1));
end
newX = s-selectedShift;
warning('on');
end