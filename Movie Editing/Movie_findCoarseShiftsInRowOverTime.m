function [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime,centerOfMass)
%[caorseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime,centerOfMass)
% Movie_findCoarseShiftsInRowOverTime will create a RowOverTime matrix
% where each row is shofted by a round number so that the asperities appear
% not to move (at least by full integers). Shifts of rows is given back by
% caorseShiftsVector.
% centerOfMass is the outpus of the function analyzeEnlargementAndStrain.

%-- create caorse shofts vector
centerOfMass = centerOfMass(:,~isnan(sum(centerOfMass)));
centerOfMassShift = centerOfMass-repmat(centerOfMass(1,:),size(centerOfMass,1),1);
centerOfMassShiftMean = mean(centerOfMassShift,2);
coarseShiftsVector = round(centerOfMassShiftMean);

%--- pad array in order to prevent data mixture from both sides of the
%matrix
maxRoundShift = max(coarseShiftsVector);
minRoundShift = min(coarseShiftsVector);
RowOverTime_padded = padarray(RowOverTime,[0 abs(minRoundShift)],0,'pre');
RowOverTime_padded = padarray(RowOverTime_padded,[0 maxRoundShift],0,'post');

%--- shoft the matrix
shiftedRowOverTime =MyCircshift(RowOverTime_padded,-coarseShiftsVector);