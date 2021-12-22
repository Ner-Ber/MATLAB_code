function bestShift = phedi_findIntrpdSlopeShiftbyMinimztn_orgnl2intrpd(slope1, slope2, numOOM, varargin)
    
    %% set parameters
    [alphaPenalty] = setDefaults4function(varargin,-1);
    
    
    %% create interpolated slopes
    simpleX = 1:length(slope1);
    intrpX = 1:10^-numOOM:length(slope1);
    try
    slope1_Intrp = interp1(simpleX,slope1,intrpX,'pchip');
    slope2_Intrp = interp1(simpleX,slope2,intrpX,'pchip');
    catch
        disp('');
    end
    weights2 = abs(gradient(slope2_Intrp)).^-alphaPenalty;
    
    %% interate on different amounts of shifts to find ideal shift
    shifts = -(2*(10^numOOM)):10^numOOM:(2*(10^numOOM));
    for OOMidx = numOOM:-1:0
        F_val = zeros(size(shifts));
        for si = 1:length(shifts)
            Idx = knnsearch([intrpX(:),slope1_Intrp(:)],[simpleX(:)+shifts(si)*10^-numOOM,slope2]);
            Di = sqrt((simpleX(:)'-intrpX(Idx)).^2 + (slope2(:)'-slope1_Intrp(Idx)).^2);
            F_val(si) = sum(weights2(Idx).*Di(:)');
        end
        [~,I] = min(F_val);
        
%         %% create infastructure for compuation with matrices
%         %     intrpX = 1:length(slope1);
%         allXs = repmat([intrpX(:); intrpX(:)],1,length(shifts));
%         shiftsMat = repmat(shifts(:)',length(slope1_Intrp),1);
%         allShiftes = [zeros(size(shiftsMat)); shiftsMat];
%         XandShifted = allXs + allShiftes;   % all x coordinates considering shifts
%         
%         %% create matrix of coordinate triplets
%         %--- create a mtrix containing triples so that the middle value is a
%         %- point and the two side values are coordinates that define the line
%         %- corosponding to it.
%         colIdx = repmat(1:length(shifts),2*length(slope1_Intrp),1);
%         [Xordered,rowIdx] = sort(XandShifted);
%         I = sub2ind(size(rowIdx),rowIdx,colIdx);
%         allYs = repmat([slope1_Intrp(:); slope2_Intrp(:)],1,length(shifts));
%         Yordered = allYs(I);
%         N = size(Xordered,1);
%         sequence = 1:(N-2);
%         TrpltsIndxsSingleShift = repmat(permute([sequence(:),sequence(:)+1,sequence(:)+2],[1 3 2]),1,length(shifts),1);
%         
%         additionMultiplyer = (0:(length(shifts)-1))*N;
%         TrpltsIndxs = bsxfun(@plus,TrpltsIndxsSingleShift,additionMultiplyer);
%         
%         allTrpletsX = Xordered(TrpltsIndxs);
%         allTrpletsY = Yordered(TrpltsIndxs);
%         u1i = allTrpletsX(:,:,1);
%         v1i = allTrpletsY(:,:,1);
%         xi = allTrpletsX(:,:,2);
%         yi = allTrpletsY(:,:,2);
%         u2i = allTrpletsX(:,:,3);
%         v2i = allTrpletsY(:,:,3);
%         
%         mi = (v2i-v1i)./(u2i-u1i);
%         Di = abs(mi.*xi-yi+v2i-u2i.*mi)./sqrt(mi.^2 + 1)./(N-2);
%         F_val = sum((abs(mi).^(-alphaPenalty)).*Di,1);
%         
%         %-- fix for full pixel shifts
%         F_valFix = F_val;
%         n = length(slope1_Intrp);
%         fullShifts = mod(shifts,1)==0;
%         fullShiftVal = shifts(fullShifts);
%         fullShiftsIdx = find(fullShifts);
%         for i=1:length(fullShiftVal)
%             F_valFix(fullShiftsIdx(i)) = rms(slope2_Intrp((max(1-fullShiftVal(i),1)):n)-slope1_Intrp(1:(min(n+fullShiftVal(i),n))));
%         end
%         
%         
%         
%         [~,I] = min(F_valFix);
        %     bestShift_i = shifts(I);
        
        %% create new shifts vector
        prev_shifts = shifts;
        bestShift_i = shifts([(I-1),(I+1)]);
        shifts = bestShift_i(1):(10^(OOMidx-1)):bestShift_i(2);
    end
    bestShift = prev_shifts(I);
end