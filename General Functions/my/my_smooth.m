function [v_s,t]=my_smooth(v,smt,t)
%The function smoothes the data including matrices, column oriented.
%if 't' is input. the edges are cutted 
%smt may be a vector

if (length(smt)==1)
    smt=ones(1,length(v(1,:)))*smt;
end

v_s=zeros(size(v));
for j=1:length(v(1,:))
    v_s(:,j)=smooth(v(:,j),smt(j));
end


if nargin>2
    v_s=v_s(smt:end-smt,:);
    t=t(smt:end-smt);
end