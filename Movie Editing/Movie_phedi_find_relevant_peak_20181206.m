function [chosenPeakLoc , side] = Movie_phedi_find_relevant_peak_20181206(enterSignal, PhediCoor)
% [chosenPeakLoc , side] = Movie_phedi_find_relevant_peak(enterSignal, PhediCoor)


        %% by peak prominance:
        %         %--- find wether there is noise
        %         NoNoise = 0;
        %         for i = 6:length(enterSignal)
        %             if isequal(diff(enterSignal((i-5):i)),zeros(size(diff(enterSignal((i-5):i)))))
        %                 NoNoise = 1;
        %                 break
        %             end
        %         end
        %
        %         [pks,locs,w,p] = findpeaks(enterSignal);
        %
        %         n = length(BOSC);
        %         m = length(locs);
        %         BOSC_duplicateMatrix = repmat(BOSC,[m 1]);
        %         locs_duplicateMatrix = repmat(locs(:),[1 n]);
        %
        %         if NoNoise
        %             locsDist = abs(locs_duplicateMatrix - BOSC_duplicateMatrix);
        %             [~,II] = min(locsDist,[],1);
        %             chosenPeakLoc = locs(II);
        %         else
        %             locsMinus = locs_duplicateMatrix - BOSC_duplicateMatrix;
        %             locsMinus(locsMinus<=0) = inf;
        %             closestFromRight= min(locsMinus,[],1)+BOSC;
        %
        %             locsMinus = BOSC_duplicateMatrix - locs_duplicateMatrix;
        %             locsMinus(locsMinus<=0) = inf;
        %             closestFromLeft = BOSC-min(locsMinus,[],1);
        %
        %             %--- get peak with higher prominence:
        %             closestPeaks = [closestFromLeft; closestFromRight];
        %             [~,Locb] = ismember(closestPeaks,locs);
        %             prominanceOfSelectedLocs = p(Locb);
        %             [BetterPromVal,~] = max(prominanceOfSelectedLocs,[],1);
        %             [~, IdxsOfBetterProm] = ismember(BetterPromVal,p);
        %             chosenPeakLoc = locs(IdxsOfBetterProm);
        %         end
        %         %--- mark negative or positive slope
        %         side = sign(chosenPeakLoc-BOSC(:));  % -1 for negative slope, +1 for positive
        
        %% by marking a slope:
        %--- find another point above this BOSC:
        side = zeros(size(PhediCoor));
        neighborFact = 1;
        while ~~nnz(~side)
            neighbors2TheRight = PhediCoor+neighborFact;
            side(side==0) = sign(enterSignal(neighbors2TheRight(side==0))-enterSignal(PhediCoor(side==0)));
            neighborFact = neighborFact+1;
        end
        
        
        [~,locs,~,~] = findpeaks(enterSignal);
        PhediCoor = PhediCoor(:)';
        locs = locs(:);
        %--- add bluff peaks at beggining and end to always ensure a results
        if min(locs)>1; locs = [1; locs]; end
        if max(locs)<length(enterSignal); locs = [locs; length(enterSignal)]; end
        
        chosenPeakLoc = zeros(size(PhediCoor));
        %--- find closest for positive slopes:
        positiveBOSC = PhediCoor(side==1);
        BOSC_mat = repmat(positiveBOSC,[length(locs), 1]);
        if ~isempty(BOSC_mat)
            locs_mat = repmat(locs,[1, length(positiveBOSC)]);
            Diff_mat = locs_mat - BOSC_mat;
            Diff_mat(Diff_mat<=0) = nan;
            chosenPeakLoc(side==1) = positiveBOSC+abs(min(Diff_mat,[],1));
        end
        
        %--- find closest for negative slopes:
        negativeBOSC = PhediCoor(side==-1);
        BOSC_mat = repmat(negativeBOSC,[length(locs), 1]);
        if ~isempty(BOSC_mat)
            locs_mat = repmat(locs,[1, length(negativeBOSC)]);
            Diff_mat = BOSC_mat - locs_mat;
            Diff_mat(Diff_mat<=0) = nan;
            chosenPeakLoc(side==-1) = negativeBOSC-abs(min(Diff_mat,[],1));
        end
    end