function [SyyOffset,Anorm]=calibrateAfromSyy(Syy,A)

A0=60;

[~,tindex]=min(abs(A(1:end,:)-A0));
for j=1:length(tindex)
    SyyOffset(:,j)=Syy(tindex(j),j);
end
 
%   Syy=Syy-repmat(SyyOffset,length(Syy(:,1)),1);
%   minSyy=min(max(Syy));

minSyy=4;
[~,tindex]=min(abs(Syy-minSyy));
for j=1:length(tindex)
    Anorm(:,j)=A(tindex(j),j);
end

% A=A./repmat(Anorm,length(A(:,1)),1);
% figure;
% plot(Syy,A,'.-');


