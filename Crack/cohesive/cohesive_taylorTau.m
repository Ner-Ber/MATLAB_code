function Y = cohesive_taylorTau(X,prefacVec)
    
    tau_power = @(x,p) (1+x).^p;
    Handle = @(x) 0;
    for i=1:length(prefacVec)
        Handle =@(x) Handle(x) + prefacVec(i)*tau_power(x,i);
    end
    HandleNorm =@(x)  Handle(x)./sum(prefacVec);
    
    Y = HandleNorm(X);
end