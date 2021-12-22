function meanDepth = profilometer_calcMeanDepthOfScratch()

Dat = data_tip_get_from_mesh;
[~, I] = sort(Dat.x);
yy = Dat.y;
yy = yy(I);
DepthDiff = diff(yy);
meanDepth = mean(abs(DepthDiff (1:2:end)));
end