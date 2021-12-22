

idx=8; % SG index

deltax=-0.015;
deltaUxx=0;
deltaUxy=-2e-4;
deltaUyy=0;
x_end=50;
GemmaVec=0.7;%0.1:0.5:2;

% [~, ~, Cr, ~, ~, ~, ~, ~, ~, ~,]=CrackSolutionMaterialProperties;
Cr=1257;
dx=1E-5;
x=-x_end*1E-3:dx:x_end*1E-3;

[~,IdxMin] = min(abs(PhediStruct.BigPicRotStruct.x-PhediStruct.SgData.x_sg(idx)/1000));
v = PhediStruct.BigPicRotStruct.frontVel_interp(IdxMin)*PhediStruct.BigPicRotStruct.fps*PhediStruct.BigPicRotStruct.res;
h = PhediStruct.SgData.y_sg(idx);
% v=GeneralStruct.Front.Vsg(idx);
% h=GeneralStruct.sg_data.y_sg(idx);

figure
subplot(1,3,1);hold all
plot(PhediStruct.SgData.x_mins_x_tips(:,idx),bsxfun(@minus,PhediStruct.SgData.Uxx(:,idx),PhediStruct.SgData.Uxx(1,idx))/1000);
for Gamma=GemmaVec;
    sol=CrackSolutionForh_GammaChange(v/Cr,-0.5,h*1e-3,Gamma,x);
    plot(sol.x-deltax,sol.Uxx+deltaUxx,'DisplayName',['\Gamma= ' num2str(sol.Gamma)])
end
title(['Uxx at SG ' num2str(idx)])
xlim([-0.1 0.1]); ylim([-0.2 0.2]/1000);

subplot(1,3,2);hold all
plot(PhediStruct.SgData.x_mins_x_tips(:,idx),bsxfun(@minus,PhediStruct.SgData.Uxy(:,idx),PhediStruct.SgData.Uxy(1,idx))/1000);
for Gamma=GemmaVec;
    sol=CrackSolutionForh_GammaChange(v/Cr,-0.5,h*1e-3,Gamma,x);
    plot(sol.x-deltax,sol.Uxy+deltaUxy,'DisplayName',['\Gamma= ' num2str(sol.Gamma)])
end
title(['Uxy at SG ' num2str(idx)])
xlim([-0.1 0.1]); ylim([-0.2 0.2]/1000);

subplot(1,3,3);hold all
plot(PhediStruct.SgData.x_mins_x_tips(:,idx),bsxfun(@minus,PhediStruct.SgData.Uyy(:,idx),PhediStruct.SgData.Uyy(1,idx))/1000);
for Gamma=GemmaVec;
    sol=CrackSolutionForh_GammaChange(v/Cr,-0.5,h*1e-3,Gamma,x);
    plot(sol.x-deltax,sol.Uyy+deltaUyy,'DisplayName',['\Gamma= ' num2str(sol.Gamma)])
end
title(['Uyy at SG ' num2str(idx)])
xlim([-0.1 0.1]); ylim([-0.2 0.2]/1000);

