function Y = fresnelComplexFun(X,varargin)
    
    [order] = setDefaults4function(varargin,55);
    
    Y = (1+1i)./sqrt(8)-(pi^-0.5).*(my_fresnelc(X,order)+1i*my_fresnels(X,order));
    
end