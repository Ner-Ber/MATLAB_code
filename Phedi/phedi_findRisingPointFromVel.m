function S = phedi_findRisingPointFromVel(Vel,Xmin_xtip,varargin)
    
    %% set defaults
    [sectionBoundary] = setDefaults4function(varargin,0.05);
    
    %% get relevant section of kink
    sectionLogic = Xmin_xtip<=sectionBoundary & Xmin_xtip>0;
    Xsection = Xmin_xtip(sectionLogic);
    Vsection = Vel(sectionLogic);
    %-- rescale Vel
    v0 = max(Vsection);
    %     Vsection = Vsection./v0;
    
    
    %% iterate on phedis to find fits
    %--- set initial and boudaries
    Eps = 1e-3;
    a0 = [sectionBoundary/10, sectionBoundary/5, v0/10];
    ub = [sectionBoundary-Eps, sectionBoundary, v0*(1-Eps)];
    lb = [Eps, 2*Eps, v0*Eps];
    %     model = @(a,x) my_Heaviside(a(1)-x).*(((a(3)-v0)./(a(1))).*x+v0) +...
    %         my_Heaviside(x-a(1)).*my_Heaviside(a(2)-x).*((x-a(2)).*a(3)./(a(1)-a(2)));
    model = @(a,x) heaviside(a(1)-x).*(((a(3)-v0)./(a(1))).*x+v0) +...
        heaviside(x-a(1)).*heaviside(a(2)-x).*((x-a(2)).*a(3)./(a(1)-a(2)));
    
    %-- exclude nans
    exNanLogic = ~(isnan(Xsection) | isnan(Vsection));
    
    opts = optimset('Display','off');
    a = lsqcurvefit(model,a0,Xsection(exNanLogic),Vsection(exNanLogic),lb,ub,opts);
    
    %% output
    S.risingVelPoint = a(2);
    S.VelKinkX = a(1);
    S.VelKinkY = a(3);
    S.v0 = v0;
    
    
end
