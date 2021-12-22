function anls_SuperShear_EqmPlot(anl)



[~, index_sg1]=min(abs(anl.y3_5.x_sg -50));
[~, index_sg2]=min(abs(anl.y3_5.x_sg -150));

index_e=[1:length(anl.e)];
anl.y3_5.Cf=zeros(length(anl.e),length(anl.y3_5.x_sg));
l=120; %location of Cf


for j=1:length(anl.e)
    %--- get Cf at sg location
    
    
        [~, index]=min(abs(anl.frontRaw{j}.xv - l));
        Cf(j)=mean( anl.frontRaw{j}.v(index-3:index+3));
        Uxy(j)=mean(anl.y3_5.Uxy0_f(j,index_sg1:index_sg2));
        
    
end

%%%%---------------------------Plot amplitudes averaged over sg indicated by index_sg


plot(4.233*Uxy,Cf,'o');
my_legend_add([anl.lgnd{1}(1:end-3) ', l=' num2str(l) ', x_sg=' num2str(anl.y3_5.x_sg(index_sg1)) '-' num2str(anl.y3_5.x_sg(index_sg2))]);
hold all;

% for j=1:length(Cf) plot(Cf(j),Uxy(j),'.'); hold all; end
% my_legend_add(num2str(anl.e'));legend off;


