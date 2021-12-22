function pairsMatrix = allPairs(A, B)
%allPairs returns all possible pairs of elements in A and B matrices. All
%inputs are treated as column vectors no matter their domensions. 
% pairsMatrix - is a (n*m)x2 matric where each row contains a pair. Here
% n,m represent the number of element in A and B.

A = A(:); B = B(:); 
C = kron(A,ones(size(B)));
D = repmat(B,length(A),1);
pairsMatrix = [C D];


end