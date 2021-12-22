function interpolate_smoothCubicSpline(signal,interpFactor,Ninterp)
    % interpolate_smoothCubicSpline(signal,interpFactor,Ninterp)
    %
    % interpolate_smoothCubicSpline is meant to create a smooth interpolant
    % by creating a moving smoothing interpolation.
    %   THIS DOESN'T WORK WELL AT ALL
    
    
    Nwindows = length(signal)-Ninterp+1;
    [Xi,YYi] = meshgrid(1:Ninterp,1:Nwindows);
    [Xq,YYq] = meshgrid(1:1/interpFactor:(Ninterp+1-1/interpFactor),1:Nwindows);
    Xinitial_mat = bsxfun(@plus,Xi,(0:(Nwindows-1))');
    Vinitial = signal(Xinitial_mat);
    partsInterped = interp2(Xi,YYi,Vinitial,Xq,YYq,'spline');
    
    Xinitial_mat_qq = bsxfun(@plus,Xq,(0:(Nwindows-1))');
    partsInterped = nan(size(Xq));
    for i = 1:size(Xinitial_mat,1)
        partsInterped(i,:) = slmeval(Xinitial_mat_qq(i,:),slmengine(Xinitial_mat(i,:),Vinitial(i,:),'degree',3,'knots',Ninterp*interpFactor));
    end
    
    
    
    NN = nan(Nwindows,length(signal)*interpFactor);
    NN(:,1:Ninterp*interpFactor) = partsInterped;
    CircShiftsVec = (0:(Nwindows-1))*interpFactor;
    NN_shifted = MyCircshift(NN,CircShiftsVec(:));
    interpedCurve = mean(NN_shifted,1,'omitnan');
    
    
end