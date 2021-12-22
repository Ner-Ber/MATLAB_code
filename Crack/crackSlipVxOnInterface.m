function [V,X] = crackSlipVxOnInterface(v,Gamma,varargin)
    
    [X] = setDefaults4function(varargin,-0.05:1e-4:0.05);
    
    [Cd, Cs, ~,~,~,~, mu, ~,~,~,~]=CrackSolutionMaterialProperties;
    G=Gamma;
    
    k=Cs/Cd;%Broberg p.330
    alpha_d=(1-(v./Cd).^2).^0.5;
    alpha_s=(1-(v./Cs).^2).^0.5;
    D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
    A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2;
    K=(G*mu*4*(1-k^2)./A).^0.5;
    
    if numel(v)==1
        V = zeros(size(X));
        V(X<=0) = (-2*pi*X(X<=0)).^-0.5.*v.*alpha_s.*K.*(1-alpha_s.^2)./mu./D;
    else
%         V = zeros(size(v));
        V = (2*pi*X*sign(X)).^-0.5.*v.*alpha_s.*K.*(1-alpha_s.^2)./mu./D;
    end
    
    
    
    
end