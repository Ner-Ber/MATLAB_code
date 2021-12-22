function bestShift = phedi_findSlopeShiftbyMinimztn(slope1, slope2, shifts, varargin)
    
    %% set parameters
    [alphaPenalty] = setDefaults4function(varargin,1e-3);
    
    %% create infa for compuation with matrices
    Xoriginal = 1:length(slope1);
    allXs = repmat([Xoriginal(:); Xoriginal(:)],1,length(shifts));
    shiftsMat = repmat(shifts(:)',length(slope1),1);
    allShiftes = [zeros(size(shiftsMat)); shiftsMat];
    XandShifted = allXs + allShiftes;   % all x coordinates considering shifts
    
    %% create matrix of coordinate triplets
    %--- create a mtrix containing triples so that the middle value is a
    %- point and the two side values are coordinates that define the line
    %- corosponding to it.
    colIdx = repmat(1:length(shifts),2*length(slope1),1);
    [Xordered,rowIdx] = sort(XandShifted);
    I = sub2ind(size(rowIdx),rowIdx,colIdx);
    allYs = repmat([slope1(:); slope2(:)],1,length(shifts));
    Yordered = allYs(I);
    N = size(Xordered,1);
    sequence = 1:(N-2);
    TrpltsIndxsSingleShift = repmat(permute([sequence(:),sequence(:)+1,sequence(:)+2],[1 3 2]),1,length(shifts),1);
    
    additionMultiplyer = (0:(length(shifts)-1))*N;
    TrpltsIndxs = bsxfun(@plus,TrpltsIndxsSingleShift,additionMultiplyer);
    
    allTrpletsX = Xordered(TrpltsIndxs);
    allTrpletsY = Yordered(TrpltsIndxs);
    u1i = allTrpletsX(:,:,1);
    v1i = allTrpletsY(:,:,1);
    xi = allTrpletsX(:,:,2);
    yi = allTrpletsY(:,:,2);
    u2i = allTrpletsX(:,:,3);
    v2i = allTrpletsY(:,:,3);

    mi = (v2i-v1i)./(u2i-u1i);
    Di = abs(mi.*xi-yi+v2i-u2i.*mi)./sqrt(mi.^2 + 1)./(N-2);
    F_val = sum((abs(mi).^(-alphaPenalty)).*Di,1); 
    
    %-- fix for full pixel shifts
    F_valFix = F_val;
    n = length(slope1);
    fullShifts = mod(shifts,1)==0;
    fullShiftVal = shifts(fullShifts);
    fullShiftsIdx = find(fullShifts);
    for i=1:length(fullShiftVal)
        F_valFix(fullShiftsIdx(i)) = rms(slope2((max(1-fullShiftVal(i),1)):n)-slope1(1:(min(n+fullShiftVal(i),n))));
    end
    
   
    
    [~,I] = min(F_valFix);
    bestShift = shifts(I);
    
end