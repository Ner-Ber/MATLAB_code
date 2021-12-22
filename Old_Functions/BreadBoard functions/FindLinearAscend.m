function ascendRange = FindLinearAscend(v, d)
%FindLinearAscend will find the linear ascending region in the calibration
%graph. The function is ment to be implemented in the "Calibration1.m"
%function.
% v - voltage vector of calibration. 1D
% d - distance vector of calibration correlating to v. 1D.


%% smooth out noise on d vector
%--- most of the noise is on the distance vector. experiment has shown that
%a good smooth is averaging upon 10% of the vector. this parameter may be
%changed.
avgDist = 0.02;
avgDistParm = round(length(d)*avgDist);
if mod(avgDistParm,2)==0; avgDistParm = avgDistParm+1; end;     % create a odd smoothing parameter
d_smooth = smooth(d,avgDistParm);






end