function [bestShift, RMS] = phedi_findIntrpdSlopeShiftbyMinimztnRMS(slope1, slope2, numOOM, varargin)
    
    %% set parameters
    [alphaPenalty] = setDefaults4function(varargin,0);
    
    
    %% create interpolated slopes
    simpleX = 1:length(slope1);
    intrpX = 1:10^-numOOM:length(slope1);
    try
    slope1_Intrp = interp1(simpleX,slope1,intrpX,'pchip');
    slope2_Intrp = interp1(simpleX,slope2,intrpX,'pchip');
%     slope1_Intrp = interp1(simpleX,slope1,intrpX,'spline');
%     slope2_Intrp = interp1(simpleX,slope2,intrpX,'spline');
    catch
        disp('');
    end
    L = length(slope1_Intrp);
    
    %% weight signals by slopes
%     F = floor(intrpX);
%     inclines1 = diff(slope1); inclines1(length(slope1)) = inclines1(length(inclines1));
    weights1 = abs(gradient(slope1_Intrp)).^-alphaPenalty;
%     inclines2 = diff(slope2); inclines2(length(slope2)) = inclines2(length(inclines2));
%     weights2 = abs(gradient(slope2_Intrp)).^-alphaPenalty;
%     slope1_Intrp_wghtd = weights1(:).*slope1_Intrp(:);
%     slope2_Intrp_wghtd = weights2(:).*slope2_Intrp(:);
    
    %% interate on different amounts of shifts to find ideal shift
    shifts = -(2*(10^numOOM)):10^numOOM:(2*(10^numOOM));
    for OOMidx = numOOM:-1:0
        %-- get indexes of points to use
        boundaries1 = [1+max(shifts(:),0),L+min(shifts(:),0)];
        boundaries2 = [1-min(shifts(:),0),L-max(shifts(:),0)];
        rmsVec = zeros(size(shifts));
        for si = 1:length(shifts)
%             rmsVec(si) = rms(slope1_Intrp_wghtd(boundaries1(si,1):boundaries1(si,2))-slope2_Intrp_wghtd(boundaries2(si,1):boundaries2(si,2)));
%             rmsVec(si) = rms(slope1_Intrp(boundaries1(si,1):boundaries1(si,2))-slope2_Intrp(boundaries2(si,1):boundaries2(si,2)));
            Nmin1 = (boundaries1(si,2)-boundaries1(si,1)+1)^-1;
            W = weights1(boundaries1(si,1):boundaries1(si,2));
            rmsVec(si) = sqrt(Nmin1.*sum(W.*abs(slope1_Intrp(boundaries1(si,1):boundaries1(si,2))-slope2_Intrp(boundaries2(si,1):boundaries2(si,2))).^2));
        end
        [minRMS,I] = min(rmsVec);
        %% create new shifts vector
        prev_shifts = shifts;
        bestShift_i = shifts([max((I-1),1),min((I+1),length(shifts))]);
        shifts = bestShift_i(1):(10^(OOMidx-1)):bestShift_i(2);
    end
    bestShift = prev_shifts(I);
    if nargout>1
        RMS = minRMS;
    end
end