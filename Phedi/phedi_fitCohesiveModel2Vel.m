function [L_chosen, Dx] = phedi_fitCohesiveModel2Vel(Cf,Gamma,Vel,x_min_x_tip,cohesiveFun,varargin)
    %% set defaults and parameters
    [verticalTolerance] = setDefaults4function(varargin,1e-4);
    [~, ~, Cr]=CrackSolutionMaterialProperties;
    Nx = 5;
    
    %% get parameters from actual data
    Xrange = [-0.06, 0.06];
    X_logical = x_min_x_tip>=Xrange(1) & x_min_x_tip<=Xrange(2);
    relevantX = x_min_x_tip(X_logical);
    relevantV = Vel(X_logical);
    [Y_peakData, I_data] = max(relevantV);
    X_peakData = relevantX(I_data);
    
    %% find parameters for fit
    DY = inf;
    L_vec = linspace(1e-6,20e-3,Nx);
    Itr = 0;
    while DY > verticalTolerance && Itr<80
        dDy = inf*ones(1,length(L_vec));
        for i = 1:length(L_vec)
            [~, Y_peak_tild] = cohesive_findPeakOfVelByIterations(Cf/Cr,1e-8,Gamma,'L',L_vec(i),cohesiveFun);
            Y_peak = -Cf*Y_peak_tild;
            dDy(i) = abs(Y_peakData-Y_peak);
        end
        [DY,I] = min(dDy);
        L_chosen = L_vec(I);
        if I==1
            L_vec = linspace(2*L_vec(1)-L_vec(2),L_vec(2),Nx);
        elseif I==length(dDy)
            L_vec = linspace(L_vec(end-1),2*L_vec(end)-L_vec(end-1),Nx);
        else
            L_vec = linspace(L_vec(I-1),L_vec(I+1),Nx);
        end
        Itr = Itr+1;
    end
    
    %% FOR
    %     L_vec = linspace(1e-6,20e-3,300);
    %     dDy = nan(1,length(L_vec));
    %     Y_peakVec = nan(1,length(L_vec));
    %     for i = 1:length(L_vec)
    %         [~, Y_peak] = cohesive_findPeakOfVelByIterations(Cf/Cr,1e-8,Gamma,'L',L_vec(i),cohesiveFun);
    %         dDy(i) = abs(Y_peakData-Y_peak);
    %         Y_peakVec(i) = Y_peak;
    %     end
    %     [DY,I] = min(dDy);
    %     L_chosen = L_vec(I);
    %     if I==1
    %         L_vec = linspace(2*L_vec(1)-L_vec(2),L_vec(2),Nx);
    %     elseif I==length(dDy)
    %         L_vec = linspace(L_vec(end-1),2*L_vec(end)-L_vec(end-1),Nx);
    %     else
    %         L_vec = linspace(L_vec(I-1),L_vec(I+1),Nx);
    %     end
    
    %%
    [X_peak, Y_peak] = cohesive_findPeakOfVelByIterations(Cf/Cr,1e-8,Gamma,'L',L_chosen,cohesiveFun);
    Dx = X_peakData - X_peak;
end