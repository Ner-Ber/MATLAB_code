function [output, lag] = my_xcorr(A, B)
%discreteXcorrelation will execute a discrete cross correlation by discrete
%convolution of two 1D input signals
%
% A,B - are the signals to check correlation upon.
% shifts2scan - limit the shifts you want to check the scan. Is optional.
%
% essentially the function will give how much A is shifted in comparison to
% B (in this order).

    [~, maxSigIdx] = max([length(A),length(B)]);
    if maxSigIdx==1
        LongSig = A;
        ShortSig = B;
        flipFlag = 1;
    elseif maxSigIdx==2
        LongSig = B;
        ShortSig = A;
        flipFlag = 0;
    end
    ShortSig = ShortSig(:)';
    LongSig = LongSig(:)';
    shortLength = length(ShortSig);
    longLength = length(LongSig);
    lag = -longLength:longLength;       % number of scans needed to X correlate
    zerosBase = zeros(1,3*longLength);    % preperation for moving signal
    output = zeros(1,length(lag));
    shortBase = [zeros(1,longLength), ShortSig, zeros(1,2*longLength-shortLength)]; % short signal with zero padding, for dot product

    %--- compute convolution
    for n = lag
        lagLocation = find(lag==n);
        LongPadded = zerosBase;
        LongPadded(lagLocation:(longLength+lagLocation-1)) = LongSig;
        output(lagLocation) = dot(LongPadded, shortBase);

    end
    if flipFlag
        output = fliplr(output);
    end

end