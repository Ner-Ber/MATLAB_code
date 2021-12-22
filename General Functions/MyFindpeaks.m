function [pks,locs,w,p] = MyFindpeaks(SignalMatrix)
%[pks,locs,w,p] = my_findpeaks(SignalMatrix)
% 
% my_findpeaks will do the same as findpeaks but for a matrix where each
% column is a signal which you would like to examine. 
% all outputs are matrixes the size of the input matrix, with NaNs in the
% places where there is now relevant info.


S = size(SignalMatrix);
N = S(2);
[pks,locs,w,p] = deal(nan(S),nan(S),nan(S),nan(S));
for i = 1:N
    [pks_i,locs_i,w_i,p_i] = findpeaks(SignalMatrix(:,i));
    pks(locs_i,i) = pks_i;
    locs(locs_i,i) = locs_i;
    w(locs_i,i) = w_i;
    p(locs_i,i) = p_i; 
end