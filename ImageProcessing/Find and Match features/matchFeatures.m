function [ind1,ind2] = matchFeatures(desc1,desc2,minScore)
% MATCHFEATURES Match feature descriptors in desc1 and desc2.
% Arguments:
% desc1 ? A kxkxn1 feature descriptor matrix.
% desc2 ? A kxkxn2 feature descriptor matrix.
% minScore ? Minimal match score between two descriptors required to be
% regarded as matching.
% Returns:
% ind1,ind2 ? These are m?entry arrays of match indices in desc1 and desc2.
%
% Note:
% 1. The descriptors of the ith match are desc1(ind1(i)) and desc2(ind2(i)).
% 2. The number of feature descriptors n1 generally differs from n2
% 3. ind1 and ind2 have the same length.


    %% evaluate score of each descriptor multiplication & creating a matrix containing scores
    % this matrix is created by matrix multiplication. each cell in the
    % resulting matrix will be the dot product of a multiplication from desc1
    % and desc2. the rows and columns indecies will indicate on the location in
    % the desc1/2 array.
    S1 = size(desc1);
    S2 = size(desc2);
    desc1_per = permute(reshape(desc1,1,S1(1)*S1(2), []),[3 2 1]);
    desc2_per = permute(reshape(desc2,1,S2(1)*S2(2), []),[2 3 1]);
    S_mat = desc1_per*desc2_per;
    S_mat(isnan(S_mat)) = -inf;
    Sdim = size(S_mat);

    %% find correct indecies according to conditions
    [~,indicesCols] = sort(S_mat,1,'descend');  % gives coordinates of largest elements in columns
    [~,indicesRows] = sort(S_mat,2,'descend');  % gives coordinates of largest elements in rows

    % [row, col] of max & 2ndmax indecies in S_mat
    max2indicesCols = [indicesCols(1,:)' , (1:Sdim(2))'; indicesCols(2,:)' , (1:Sdim(2))']; 
    max2indicesRows = [(1:Sdim(1))' , indicesRows(:,1); (1:Sdim(1))', indicesRows(:,2) ]; 

    max2indicesCols_i = sub2ind(Sdim,max2indicesCols(:,1),max2indicesCols(:,2));
    max2indicesRows_i = sub2ind(Sdim,max2indicesRows(:,1),max2indicesRows(:,2));

    BothMax = max2indicesCols_i(ismember(max2indicesCols_i, max2indicesRows_i));
    AboveThresh = find(S_mat >= minScore);

    [ind1,ind2] = ind2sub(Sdim,BothMax(ismember(BothMax,AboveThresh))); 

end