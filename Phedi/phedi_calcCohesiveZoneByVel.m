function x_c = phedi_calcCohesiveZoneByVel(PhediVelocity,x_mins_x_tip)
% x_c = phedi_calcCohesiveZoneByVel(PhediVelocity,x_mins_x_tip)
%
% phedi_calcCohesiveZoneByVel - will calculate the cihesive zone by
% assuming that the velocity v(x)=C/sqrt(x) on the interface. By comparing
% the v(x_0), x_c (cohsive length) can be determined.
% this will be done by the calculation:
% x_c=epsilon/((v_Xc/v_Xc_ans_epsln)^2-1).
%
% PhediVelocity - a matrix containing phedi velocity vectors. Each column
% is the velocity for a different phedi. 
% x_mins_x_tip - location vector where the location of the front is at x=0.

examineLength = 1.5*1e-3; % the length (epsilon) upon which the reference calculation will be made (then averaged)
%--- search in region
searchLength = 1.5e-3;
[~,beginBound] = min(abs(x_mins_x_tip+searchLength));
[~,endBound] = min(abs(x_mins_x_tip-searchLength));

[v_Xc,I_Xc] = max(PhediVelocity(beginBound:endBound,:),[],1);
I_Xc = I_Xc+beginBound-1;
x_c = nan(1,size(PhediVelocity,2));
for i = 1:size(PhediVelocity,2)
    x_relative2peakVel = x_mins_x_tip - x_mins_x_tip(I_Xc(i));
    relevantXindexes = x_relative2peakVel>0 & x_relative2peakVel<=examineLength;
    epsilon = x_relative2peakVel(relevantXindexes);
    x_c_vec = epsilon./((v_Xc(i)./PhediVelocity(relevantXindexes(:),i)).^2-1);
    x_c(i) = mean(x_c_vec);
end

end