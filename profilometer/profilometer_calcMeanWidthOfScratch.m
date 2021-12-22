function meanWidth = profilometer_calcMeanWidthOfScratch()

Dat = data_tip_get_from_mesh;
XDat = sort(Dat.x);
WidthDiff = diff(XDat);
meanWidth = mean(WidthDiff(1:2:end));
end