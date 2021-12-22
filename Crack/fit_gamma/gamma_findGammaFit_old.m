function [] = gamma_findGammaFit(SgDataStruct,sg_indexes)
%gamma_findGammaFit will find the gamma of a sg based on residual
%minimization

%% set values to iterate upon
% [~, ~, ~, ~, ~, ~, ~, gamma_0, ~, ~,]=CrackSolutionMaterialProperties;
Cr = 1257;
% gamma_tryVec = gamma_0-0.5:0.05:gamma_0+0.5;
gamma_tryVec  = 0.1:0.3:4;
dx_tryVec = -10:2:10;
%--- Uij are in milistrain
dUxx_tryVec = -0.5:0.07:0.5;
dUxy_tryVec = -0.5:0.07:0.5;
dUyy_tryVec = -0.5:0.07:0.5;

dx=1e-5;
x_end=50e-3;
% x=-x_end:dx:x_end;

% [~,IdxMin] = min(abs(PhediStruct.BigPicRotStruct.x-PhediStruct.SgData.x_sg(sg_indexes)/1000));
% v = PhediStruct.BigPicRotStruct.frontVel_interp(IdxMin)*PhediStruct.BigPicRotStruct.fps*PhediStruct.BigPicRotStruct.res;
v = SgDataStruct.v_sg(sg_indexes);
h = SgDataStruct.y_sg(sg_indexes);
fitMat = nan(length(gamma_tryVec),length(dx_tryVec),length(dUxx_tryVec),length(dUxy_tryVec),length(dUyy_tryVec));
relevantXX = SgDataStruct.x_mins_x_tips(:,sg_indexes);
x_logical = relevantXX<x_end & relevantXX>-x_end;

relevant_sg_Uxx = SgDataStruct.Uxx(x_logical,sg_indexes)-SgDataStruct.Uxx(1,sg_indexes);
relevant_sg_Uxy = SgDataStruct.Uxy(x_logical,sg_indexes)-SgDataStruct.Uxy(1,sg_indexes);
relevant_sg_Uyy = SgDataStruct.Uxx(x_logical,sg_indexes)-SgDataStruct.Uyy(1,sg_indexes);

for i_1 = 1:length(gamma_tryVec)
    for i_2 = 1:length(dx_tryVec)
        for i_3 = 1:length(dUxx_tryVec)
            for i_4 = 1:length(dUxy_tryVec)
                for i_5 = 1:length(dUyy_tryVec)
                    %-- define x vector to find
                    newX = relevantXX;
                    newX = newX(x_logical);
                    newX = newX+dx_tryVec(i_2);
                    %-- get relevant fucntion
                    sol=CrackSolutionForh_GammaChange(v/Cr,-0.5,h*1e-3,gamma_tryVec(i_1),newX);
                    %-- fix strains for fitting
                    newUxx = sol.Uxx+dUxx_tryVec(i_3);
                    newUxy = sol.Uxy+dUxy_tryVec(i_4);
                    newUyy = sol.Uyy+dUyy_tryVec(i_5);
                    %--- find fitting residuals Uxx
                    res_xx = sum(abs(relevant_sg_Uxx-newUxx));
                    res_xy = sum(abs(relevant_sg_Uxy-newUxy));
                    res_yy = sum(abs(relevant_sg_Uyy-newUyy));
                    fitMat(i_1,i_2,i_3,i_4,i_5) = res_xx+res_xy+res_yy;
                end
            end
        end
    end
end


