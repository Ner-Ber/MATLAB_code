function gamma_plotFrontWGamma(DataStruct)

%-- plot ROT
figure; IDT_PlotRowOverTime(DataStruct.BigPicRotStruct);
%--- determine axis
x = DataStruct.BigPicRotStruct.x;
xLogic = x>0.02 & x<0.13;
frontTime_s = DataStruct.BigPicRotStruct.frontTime_interp/DataStruct.BigPicRotStruct.fps;
MAX = max(frontTime_s(xLogic));
MIN = min(frontTime_s(xLogic));
plusDy = (MAX-MIN)/10;
ylim([MIN-plusDy MAX+plusDy]);

[X,GammaVec] = gamma_calcGammaAlongInterface(DataStruct);
hold on;
yyaxis right;
plot(X*1e-3,GammaVec,'s-','MarkerFaceColor','k');