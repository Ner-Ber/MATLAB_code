% [UxxAmp,UyyAmp,UxyAmp,Cf_vec] = CrackStrainAmpVsCf(Gamma,Cf_vec,x)

function [UxxAmp,UyyAmp,UxyAmp,Cf_vec] = CrackStrainAmpVsCf(Gamma,varargin)
    [Cf_vec,x] = setDefaults4function(varargin,1:15:1255,-0.1:5e-4:0.1);
    [Cd, Cs, Cr, nu, ro, E, mu, ~,~,~,~]=CrackSolutionMaterialProperties;
    h=3.5e-3;
    
    UxxAmp = nan(size(Cf_vec));
    UyyAmp = nan(size(Cf_vec));
    UxyAmp = nan(size(Cf_vec));
    for cf_i = 1:length(Cf_vec)
        sol = CrackSolutionForh_GammaChange(Cf_vec(cf_i)/Cr,-0.5,h,Gamma,x);
        UxxAmp(cf_i) = min(sol.Uxx);
        UyyAmp(cf_i) = max(sol.Uyy);
        UxyAmp(cf_i) = max(sol.Uxy);
    end
    
    
end