function s=structCutTime(sIn,tStart,tEnd)
% sIn-time dependent structure. like acqE,A,phE. 
%cut the time dpendent dimention . dim=1;

names=fieldnames(sIn);

[~,index1]=min(abs(sIn.t-tStart));
[~,index2]=min(abs(sIn.t-tEnd));


for j=1:length(names)
   if (size(sIn.(names{j}),1)==size(sIn.t,1))
       s.(names{j})=sIn.(names{j})(index1:index2,:);
   else
       s.(names{j})=sIn.(names{j});
    end
   
end


