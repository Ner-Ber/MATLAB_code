function [edges] = IDT_findEdgePixlesOfBlock(RowOverTime,varargin)
%edges = IDT_findEdgePixlesOfBlock(RowOverTime,thresh)
% will fins the pixles where the edges of the block are. This function was
% biuld specifically for the big-picture with the IDT camers with the lens:
% nikon 35mm (but may work in more configurations!)
%
% varargin:
% thresh default = 9e-5
% 
%
%
%   09-12-2018
%       *) the edges of the block will be defined by the first and last
%       pixles which amplitudes are larger than the STD of the noise found
%       in 'AmpRow'. 


% [thresh,findEdgeDirt] = setDefaults4function(varargin,9e-5, 20);
% [thresh,findEdgeDirt] = setDefaults4function(varargin,0.1, 20);
[thresh] = setDefaults4function(varargin,0.1);


summingRow = sum(RowOverTime);
sumMean = movmean(summingRow,10);
AmpRow = summingRow-sumMean;
STD = std(AmpRow);
smallFluctn = abs(AmpRow)>=STD;

edges = [find(smallFluctn,1,'first') find(smallFluctn,1,'last')];
% edges = [find(AmpRow>thresh,1,'first') find(AmpRow>thresh,1,'last')];


% derivSum = diff(summingRow);
% L = length(summingRow);
% [~, leftBound] = max(derivSum(1:round(0.25*L)));
% [~, rightBound] = min(derivSum(round(0.75*L):end));
% rightBound = rightBound+round(0.7*L);
% edges = [leftBound rightBound];

% ROT_norm = bsxfun(@rdivide,RowOverTime(1:100,:),RowOverTime(1,:));
% ROT_norm_mean = mean(ROT_norm(1:100,:),1);
% ROT_norm_diff = diff(ROT_norm_mean);
% ROT_norm_logical = ROT_norm_diff>thresh;
% 
% %--- left edge:
% Ledge = find(ROT_norm_logical(1:round(length(ROT_norm_logical)/3)),1,'last');
% if isempty(Ledge)
%     Ledge = 1;
% end
% %--- right edge:
% Redge = length(ROT_norm_logical)-find(fliplr(ROT_norm_logical(round(2*length(ROT_norm_logical)/3):end)),1,'last')+1;
% if isempty(Redge)
%     Redge = length(ROT_norm_logical);
% end
% 
% edges = [Ledge Redge];


% 
% kernel = [-1 0 1];
% EdgesMap = conv2(RowOverTime,kernel,'same');
% EdgesMap = EdgesMap>thresh;
% 
% Ed_sum = sum(EdgesMap,1);
% Ed_sum(1) = 0;
% Ed_sum(end) = 0;
% 
% edges = [find(Ed_sum>findEdgeDirt,1,'first') find(Ed_sum>findEdgeDirt,1,'last')];
 

end

