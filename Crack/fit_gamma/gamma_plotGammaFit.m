function gamma_plotGammaFit(DataStruct,varargin)


NsgDef = DataStruct.solAtSG.SG_calc;
Nsg = setDefaults4function(varargin,NsgDef);
figure; hold on;
%-- plot nans for legend
plot(nan,nan,'Color',rgb('Red'),'LineWidth',1.5);
plot(nan,nan,'Color',rgb('Gold'),'LineWidth',1.5);
plot(nan,nan,'Color',rgb('SteelBlue'),'LineWidth',1.5);

%--- plot Uxx
plot(DataStruct.SgData.x_mins_x_tips(:,Nsg),1e-3*(DataStruct.SgData.Uxx(:,Nsg)-mean(DataStruct.SgData.Uxx(1:100,Nsg))),...
    '.-','Color',rgb('HotPink'));
plot(DataStruct.solAtSG.x,DataStruct.solAtSG.Uxx,...
    'Color',rgb('Red'));
%--- plot Uyy
plot(DataStruct.SgData.x_mins_x_tips(:,Nsg),1e-3*(DataStruct.SgData.Uyy(:,Nsg)-mean(DataStruct.SgData.Uyy(1:100,Nsg))),...
    '.-','Color',rgb('Gold'));
plot(DataStruct.solAtSG.x,DataStruct.solAtSG.Uyy,...
    'Color',rgb('Sienna'));
%--- plot Uxy
zeroRegion = [-0.05 -0.065];
referenceSG_logical = (DataStruct.SgData.x_mins_x_tips(:,Nsg)<=zeroRegion(1)) & (DataStruct.SgData.x_mins_x_tips(:,Nsg)>=zeroRegion(2));
sgStrain = DataStruct.SgData.Uxy(:,Nsg)-mean(DataStruct.SgData.Uxy(referenceSG_logical,Nsg));
referenceSOL_logical = (DataStruct.solAtSG.x<=zeroRegion(1)) & (DataStruct.solAtSG.x>=zeroRegion(2));
solStrain = DataStruct.solAtSG.Uxy-mean(DataStruct.solAtSG.Uxy(referenceSOL_logical));
plot(DataStruct.SgData.x_mins_x_tips(:,Nsg),1e-3*sgStrain,...
    '.-','Color',rgb('Turquoise'));
plot(DataStruct.solAtSG.x,solStrain,'Color',rgb('SteelBlue'));


legend('Uxx','Uyy','Uxy');
xlabel('x-x_{tip} [m]');
ylabel('Uij [Strain]');
xlim([-0.05 0.05]);
title(['fit of angular strain function with \Gamma=',num2str(DataStruct.solAtSG.Gamma)]);
end