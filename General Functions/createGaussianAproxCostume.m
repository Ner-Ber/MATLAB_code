function g = createGaussianAproxCostume(kernelHeight,kernelLength)
    %createGaussianAprox creates a 2D binomial coefficient matrix (gaussian
    %approximation) using convolution method.
    %matrix's size will be kernelSize*kernelSize
    
    if nargin<2
        kernelLength = kernelHeight;
    end
    
    nL = kernelLength - 1;
    kL = 1:nL;
    BinomCoL = [1 cumprod((nL-kL+1)./kL)];
    
    nH = kernelHeight - 1;
    kH = 1:nH;
    BinomCoH = [1 cumprod((nH-kH+1)./kH)];
    
    g_notNorm = conv2(BinomCoL, BinomCoH');       % creating a kernel matrix
    g = g_notNorm/sum(sum(g_notNorm));          % normalizing the kerna
    
    
end