function [totX,totT] = IDT_AddAnfFixPointsOfFront(NormalizedRowOverTime,X,T,varargin)
    %new_frontStructure = IDT_AddAnfFixPointsOfFront(NormalizedRowOverTime,frontStructure,varargin)
    %
    %[smthL,minFrmDist,searchDist,FallingPercnt] = setDefaults4function(varargin,5,40,25,5);
    
    [smthL,minFrmDist,searchDist,FallingPercnt] = setDefaults4function(varargin,5,40,25,5);
    smthL = smthL+~mod(smthL,2);
    
%     X = frontStructure.StepsPix;
%     T = frontStructure.StepTime;
    
    %% improve density of points
    DT = diff(T);
    longTime = find(DT>=minFrmDist);
    newpointsX = [];
    newpointsT = [];
    for i=1:length(longTime)
        newT = T(longTime(i)):minFrmDist/2:T(longTime(i)+1);
        newT = newT(2:end-1);
        newX = round(interp1([T(longTime(i)),T(longTime(i)+1)],[X(longTime(i)),X(longTime(i)+1)],newT,'linear'));
        [newX,ia,~] = unique(newX);
        newT = newT(ia);
        X4Search = bsxfun(@plus,newX(:),(-searchDist:searchDist));
        X4Search(X4Search<1) = 1;
        X4Search(X4Search>size(NormalizedRowOverTime,2)) = size(NormalizedRowOverTime,2);
        T4search = repmat(newT(:),1,size(X4Search,2));
        linearInd = sub2ind(size(NormalizedRowOverTime), T4search, X4Search);
        parts4Search = NormalizedRowOverTime(linearInd);
        parts4SearchPad = padarray(parts4Search,[0 floor(smthL/2)],'replicate');
        parts4SearchSmth = conv2(parts4SearchPad,ones(1,smthL)/smthL,'valid');
        for j=1:size(parts4SearchSmth,1)
            [Xinter,~] = intersections(X4Search(j,:),parts4SearchSmth(j,:),[-1 1]*inf,[1 1]*(1-FallingPercnt*1e-2));
            [~,I] = min(abs(Xinter-(searchDist+1)));
            PointX = round(Xinter(I));
            if ~isempty(PointX)
                newpointsX = cat(1,newpointsX,PointX);
                newpointsT = cat(1,newpointsT,newT(j));
            end
        end
    end
    
    %-- eliminate duplicates, prefer old points:
    [~,ia,~] = intersect(newpointsX(:),X(:));
    newpointsX(ia) = [];
    newpointsT(ia) = [];
    
    [totX,I] = sort([newpointsX(:);X(:)]);
    totT_tag = [newpointsT(:);T(:)];
    totT = totT_tag(I);
    
    [totX,ia] = unique(totX);
    totT = totT(ia);
    
%     new_frontStructure = frontStructure;
%     new_frontStructure.StepsPix = totX;
%     new_frontStructure.StepTime = totT;
%     new_frontStructure.StepsPix_original = X;
%     new_frontStructure.StepTime_original = T;
    
end