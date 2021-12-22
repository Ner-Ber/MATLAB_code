function A = CrackContractionFunc(v,varargin)
    
    [PlaneFlag] = setDefaults4function(varargin,'PlaneStrain');
    [Cd, Cs, Cr, nu , ~, ~, mu, GammaDefault, PlaneStrain, tau_p_def]=CrackSolutionMaterialProperties(PlaneFlag);
    
    k=Cs/Cd;%Broberg p.330
    v=v*Cr;
    alpha_d=(1-(v/Cd).^2).^0.5;
    alpha_s=(1-(v/Cs).^2).^0.5;
    D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2; % the Rayleight function
    
    if ~strcmpi(PlaneFlag,'PlaneStrain')
        A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2; % this is plain stress %Broberg p.334,336, equation 6.2.26
    else
        A=alpha_s.*v.^2./D/(1-nu)/Cs^2; %Only for plain strain
    end
    
end