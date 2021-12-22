function g = createGaussianAprox(kernelSize)
%createGaussianAprox creates a 2D binomial coefficient matrix (gaussian
%approximation) using convolution method.
%matrix's size will be kernelSize*kernelSize

    n = kernelSize - 1;
    k = 1:n;
    BinomCo = [1 cumprod((n-k+1)./k)];
    g_notNorm = conv2(BinomCo, BinomCo');       % creating a kernel matrix
    g = g_notNorm/sum(sum(g_notNorm));          % normalizing the kerna


end