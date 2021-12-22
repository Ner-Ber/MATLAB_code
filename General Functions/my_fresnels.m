function Y = my_fresnels(X,order)
    if nargin<2
        order=10;
    end
    
    nthTerm = @(x,n) ((-1)^n*x.^(4*n+3))./(factorial(2*n+1).*(4*n+3));
    Y = 0;
    for i=0:order
        Y = Y + nthTerm(X,i);
    end
    