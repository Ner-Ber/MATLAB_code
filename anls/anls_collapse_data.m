function anls_collapse_data(t,y,frontV,t0,y0)

if nargin<5
    y0=1;
end
    
for j=1:length(t0) 
    [~,index(j)]=min(abs(t-t0(j)));
    y_norm(:,j)=circshift(y(:,j),-index(j)+index(1))./y0(j);
    t_norm(:,j)=-frontV(j)*(t-t0(j));
end
figure;
plot(t_norm,y_norm,'.-');
