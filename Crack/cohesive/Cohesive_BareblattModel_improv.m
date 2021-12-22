function Ux = Cohesive_BareblattModel_improv(x,sqrtSol,d,Amp)
    %Ux = Cohesive_BareblattModel_improv(x,sqrtSol,d,Amp)
    
    l = 1/d;
    Ux = (1-exp(l*x)).*sqrtSol + exp(l*x).*Amp.*(-x).^1.5;
    Ux(x>=0) = 0;
    
end