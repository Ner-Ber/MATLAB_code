function [PhediVelocity, timeVec_4vel] = phedi_calcVelocity(PhediLocation,timeVec)
% [PhediVelocity, timeVec_4vel] = phedi_calcVelocity(PhediLocation,timeVec)
% phedi_calcVelocity - will calculate the velocity of the phedis given in
% the matrix PhediLocation.
% PhediLocation - is found in the structure in analyzePhediCell. each
%       column represtnts the phedi location in a specific time. 
% timeVec - time corresponding to the phedi location.
%       length(timeVec)=size(PhediLocation,1) should hold.


PhediVelocity = bsxfun(@rdivide, diff(PhediLocation,1,1),diff(timeVec(:)));
timeVec_4vel = movmean(timeVec,2,'Endpoints','discard');

end