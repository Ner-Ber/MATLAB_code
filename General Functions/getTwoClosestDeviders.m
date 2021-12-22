function [m,n] = getTwoClosestDeviders(X)

r = sqrt(X);
q = floor(r);
while mod(X,q)~=0
    q = q-1;
end

m = max(q,X/q);
n = min(q,X/q);