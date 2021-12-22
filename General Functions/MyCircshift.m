function Y = MyCircshift(A,K,dim)
%Y = MyCircshift(A,K,dim)
% MyCircshiftis the same as circshift but is can use different circ values
% for each row.
% K isvector in the length of the wanted dimension to shift.
% for example, if A = [1 2 3; 1 2 3; 1 2 3; 1 2 3] and K= [1 2 3 4] so
% Y= [3 1 2; 2 3 1; 1 2 3; 3 1 2]
% K is in the length of the dimensio npependicular to the one sifthed
% around;
% dim - is the dimension shifted around. in the above example dim=2.
% default is dim=2, meaning rows are circulated.
% !!works for 2D matrixs only!!

if nargin<3
    dim=2;
end
dimRotationPerp = ~(dim-1)+1;   
N = size(A,dim);
M = size(A,dimRotationPerp);
if M~=length(K)
    error('number of elements in K doesn''t fit length of relevant dimension')
end

if dim==1
    A = A';
end

Y = A;
for i = 1:M
    relevantRow = A(i,:);
    Y(i,:) = circshift(relevantRow,[0 K(i)]);
end

if dim==1
    Y = Y';
end

end

