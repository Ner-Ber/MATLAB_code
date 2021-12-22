function [BestShift, shifted_xAxis] = phedi_findShiftByFit2PowerLaw(x,phediLoc,power,amplitude,xRange)
% phedi_findShiftByFit2PowerLaw(x,pehdiLoc,power,amplitude,xRange)
%
% phedi_findShiftByFit2PowerLaw will findthe shift in ehich the phedi
% trjectory most fits the power law of form amplitude*x^power.
% a spcific range is serched between [min(xRange) max(xRange)].


% figure; hold on;

Shifts = linspace(-0.007,0.007,100);
linearModel = @(log_x) log10(amplitude)+power.*log_x;

RMS_vec = nan(1,length(Shifts));
% COLORS = MyVaryColor(length(Shifts));
for s = 1:length(Shifts)
    ShiftedX = (x+Shifts(s));  
    
    cropLogical = ShiftedX>=min(xRange) & ShiftedX<=max(xRange);
    x_cropped = ShiftedX(cropLogical);
    phediLocCropped = phediLoc(cropLogical);
    
    x_reverse = -x_cropped; % '-' sign to bring power to positive side
    
    positiveLogical = x_reverse>0;
    positiveX = x_reverse(positiveLogical); % this is the small portion to check. shifted, and reversed to positive
    logX = log10(positiveX);
    logPhediLoc = log10(phediLocCropped(positiveLogical));
    logModel = linearModel(logX);
    
    RMS_vec(s) = rms(logPhediLoc-logModel);
    
    %     plot(logPhediLoc,logModel,'Color',COLORS(s,:));
end
% axis equal
% axis square

SMTH_RMS = supsmu(Shifts,RMS_vec);
[~,I] = min(SMTH_RMS);
BestShift = Shifts(I);
shifted_xAxis = x+BestShift;

end