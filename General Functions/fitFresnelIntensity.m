function Y = fitFresnelIntensity(x0,X)
    
    Y = x0(1).*abs(fresnelComplexFun(x0(2)*(X-x0(3)))).^2+x0(4);

end