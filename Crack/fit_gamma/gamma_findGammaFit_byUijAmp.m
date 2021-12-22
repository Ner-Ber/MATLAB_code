function minimizingGamma = gamma_findGammaFit_byUijAmp(Cf_vec,UijAmp,varargin)
    [fitCmponnt,Gamma_vec] = setDefaults4function(varargin,'Uxx',linspace(0.8,11,60));
    
    RMS_vec = inf(size(Gamma_vec));
    for g = 1:length(Gamma_vec)
        [UxxAmp,UyyAmp,UxyAmp,~] = CrackStrainAmpVsCf(Gamma_vec(g),Cf_vec);
        if strcmpi(fitCmponnt,'Uxx')
            AmpTheory = UxxAmp;
        elseif strcmpi(fitCmponnt,'Uxy')
            AmpTheory = UxyAmp;
        elseif strcmpi(fitCmponnt,'Uyy')
            AmpTheory = UyyAmp;
        end
        RMS_vec(g) = rms(AmpTheory-UijAmp);
    end
    
    [~,I] = min(RMS_vec);
    minimizingGamma = Gamma_vec(I);
end