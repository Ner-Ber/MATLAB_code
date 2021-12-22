function [chosenPeakLoc , side] = Movie_phedi_find_relevant_peak(enterSignal, PhediCoor,varargin)
% [chosenPeakLoc , side] = Movie_phedi_find_relevant_peak(enterSignal, PhediCoor)

        [smt4SideFind] = setDefaults4function(varargin,3);
        
        %% by marking a slope:
        enterSignalSmt = smooth(enterSignal,smt4SideFind);
        %--- find another point above this BOSC:
        sideChk1 = zeros(size(PhediCoor));
        neighborFactVec1 = zeros(size(PhediCoor));
        neighborFact = 1;
        while ~~nnz(~sideChk1)
            neighbors2TheRight = PhediCoor+neighborFact;
            neighborFactVec1(sideChk1==0) = neighborFact;
            sideChk1(sideChk1==0) = sign(enterSignalSmt(neighbors2TheRight(sideChk1==0))-enterSignalSmt(PhediCoor(sideChk1==0)));
            neighborFact = neighborFact+1;
        end
        
        sideChk2 = zeros(size(PhediCoor));
        neighborFactVec2 = zeros(size(PhediCoor));
        neighborFact = 1;
        while ~~nnz(~sideChk2)
            neighbors2TheRight = PhediCoor-neighborFact;
            neighborFactVec2(sideChk2==0) = neighborFact;
            sideChk2(sideChk2==0) = -sign(enterSignalSmt(neighbors2TheRight(sideChk2==0))-enterSignalSmt(PhediCoor(sideChk2==0)));
            neighborFact = neighborFact+1;
        end
        
        %-- choose the closest slope:
        [~,I] = min([neighborFactVec1(:), neighborFactVec2(:)],[],2);
        side = zeros(size(PhediCoor));
        side(I==1) = sideChk1(I==1);
        side(I==2) = sideChk2(I==2);
        
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