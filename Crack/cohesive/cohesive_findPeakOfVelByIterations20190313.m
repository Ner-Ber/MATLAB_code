function [X_peak, Y_peak] = cohesive_findPeakOfVelByIterations(v,h,Gamma,variableName1,variableVal1,varargin)
    % cohesive_findPeakOfVelByIterations(v,h,Gamma,variableName1,variableVal1,varargin)
    %
    % cohesive_findPeakOfVelByIterations will find the peak in Uxx of the
    % solution of cohesive zone model (which is the same as the peak in
    % location in particle velocity) by iterating on a low resolution solution.
    %
    % OPTIONAL = defaults:
    % [cohesiveFun,PlaneFlag,convergeRes] = ['exp','PlaneStrain',5e-6]
    %
    % the function operates on the minima of Uxx since its the same location as
    % the maxima of velocity
    
    warning('off');
    
    [cohesiveFun,PlaneFlag,convergeRes] = setDefaults4function(varargin,'exp','PlaneStrain',5e-6);
    
    Nx = 30;
    X_vec = linspace(-2,0,Nx);
    DX = X_vec(2)-X_vec(1);
    while DX>=convergeRes
        DX = X_vec(2)-X_vec(1);
        currentCohesiveStruct = CrackSolutionGeneralCohesive_InsertVariables_Neri(v,h,Gamma,variableName1,variableVal1,...
            cohesiveFun,PlaneFlag,X_vec,DX);
        [pks,locs] = findpeaks(-currentCohesiveStruct.Uxx);
        %-- if there are no peaks search in broader region
        if isempty(pks)
            %         X_vec = linspace(min(X_vec)+1,max(X_vec),Nx);
            X_vec = linspace(X_vec(2),max(X_vec),Nx);
            continue
        end
        %-- take only peaks in negative Y axis.
        negativePksLogic = pks>0;
        newRegion = sort(currentCohesiveStruct.x(([-1 +1]+locs(negativePksLogic))));
        X_vec = linspace(newRegion(1),newRegion(2),Nx);
    end
    X_peak = currentCohesiveStruct.x(negativePksLogic);
    Y_peak = -pks(negativePksLogic);
    
    warning('on');
end