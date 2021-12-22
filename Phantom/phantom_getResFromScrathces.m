function res = phantom_getResFromScrathces(RowOverTime,varargin)


%% set defaults
[scratchDist, smth, UsrCutOff] = setDefaults4function(varargin,200e-6,4,nan);

%% calc
%--- get typical row
meanRow = mean(RowOverTime(1:100,:),1);
%--- smooth and flip to find scratches
smoothedRow = smooth(meanRow,smth);
% upDwnRow = -smoothedRow;
upDwnRow = smoothedRow;
%--- find scratches by using find peaks
[~,locs,~,p] = findpeaks(upDwnRow);
% %--- choose only center third where intensity is high
% L = length(upDwnRow);
% borders = round([L/3, 2*L/3]);
% relevantPeak = (locs>=borders(1) & locs<=borders(2));
%--- choose all peaks:
relevantPeak = logical(ones(size(p)));
%--- take prominance of relevant peaks (actually scratches)
relevantP = p(relevantPeak);
relevantLocs = locs(relevantPeak);
%--- create cutoff for actual peaks (this assumes that there is noise!!)
if isnan(UsrCutOff)
%     CutOff = mean([min(relevantP), max(relevantP)]);
    CutOff = max(smoothedRow)/4;
else
    CutOff = UsrCutOff;
end
%--- take only actual peaks not noise
actualLogical = relevantP>CutOff;
actualLocs = relevantLocs(actualLogical);
%--- GET DISTANCE BETWEEN SCRATHCES
scratchDistPix = mean(diff(actualLocs));
%--- calc resolution
res = scratchDistPix/scratchDist;
end