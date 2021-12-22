function anls_SuperShearAmpPlot(anl)



index_sg=9:12; %Choose sg for plot
index_e=[1:length(anl.e)];
anl.y3_5.Cf=zeros(length(anl.e),length(anl.y3_5.x_sg));

for j=1:length(anl.e)
    %--- get Cf at sg location
    
    for l=1:length(index_sg)
        
        [~, index]=min(abs(anl.frontRaw{j}.xv - anl.y3_5.x_sg(index_sg(l))));
        if (index<4)
            index=4;
        end
        anl.y3_5.Cf(j,index_sg(l))=mean( anl.frontRaw{j}.v(index-3:index+3) );
    end
end

%%%%---------------------------Plot amplitudes averaged over sg indicated
%%%%by index_sg


%figure;
subplot(2,2,1);hold all;
% for j=1:length(anl.Cf) plot(anl.Cf(j)+0*anl.y3_5.Uxy0_f(j,:),anl.y3_5.Uxy0_f(j,:),'r.');hold all; end
% my_legend_add(anl.lgnd);
% for j=1:length(anl.Cf) plot(anl.Cf(j)+0*anl.y7_5.Uxy0_f(j,:),anl.y7_5.Uxy0_f(j,:),'b.');hold all; end
% my_legend_add(anl.lgnd);
plot(mean(anl.y3_5.Cf(index_e,index_sg),2),mean(anl.y3_5.Uxy0_f(index_e,index_sg),2),'ro');
my_legend_add(['y=3.5, x=' num2str(anl.y3_5.x_sg(index_sg))]);legend off;
% plot(anl.Cf,mean(anl.y7_5.Uxy0_f,2),'bo');
% my_legend_add(['y=7.5, x=' num2str(anl.y7_5.x_sg)]);legend off;
ylabel('Uxy');
title(anl.lgnd{1}(1:19));

subplot(2,2,2);hold all;
% for j=1:length(anl.Cf) plot(anl.Cf(j)+0*anl.y3_5.UyyP(j,:),anl.y3_5.UxxP(j,:)-anl.y3_5.Uxx1(j,:),'r.');hold all; end
% my_legend_add(anl.lgnd);
% for j=1:length(anl.Cf) plot(anl.Cf(j)+0*anl.y7_5.UyyP(j,:),anl.y7_5.UyyP(j,:)-anl.y7_5.UyyM(j,:),'b.');hold all; end
% my_legend_add(anl.lgnd);
plot(mean(anl.y3_5.Cf(index_e,index_sg),2),mean(anl.y3_5.UxxP(index_e,index_sg)-anl.y3_5.Uxx1(index_e,index_sg),2),'ro');
my_legend_add('Uxx');
plot(mean(anl.y3_5.Cf(index_e,index_sg),2),mean(anl.y3_5.UyyP(index_e,index_sg)-anl.y3_5.Uyy1(index_e,index_sg),2),'r.');
my_legend_add('Uyy');
% plot(anl.Cf,mean(anl.y7_5.UxxP-anl.y7_5.Uxx1,2),'bo');
% my_legend_add('Uxx');
% plot(anl.Cf,mean(anl.y7_5.UyyP-anl.y7_5.Uyy1,2),'b.');
% my_legend_add('Uyy');
ylabel('Us');

subplot(2,2,3);hold all;
% for j=1:length(anl.Cf) plot(anl.Cf(j)+0*anl.y3_5.UxxP(j,:),anl.y3_5.UxxP(j,:),'r.');hold all; end
% my_legend_add(anl.lgnd);
% for j=1:length(anl.Cf) plot(anl.Cf(j)+0*anl.y7_5.UxxP(j,:),anl.y7_5.UxxP(j,:),'b.');hold all; end
% my_legend_add(anl.lgnd);
plot(mean(anl.y3_5.Cf(index_e,index_sg),2),mean(anl.y3_5.UxxP(index_e,index_sg),2),'ro');
%plot(mean(anl.y7_5.Cf(index_e,index_sg),2),mean(anl.y7_5.UxxP(index_e,index_sg),2),'bo');
ylabel('UxxP');
 
subplot(2,2,4);hold all;
% for j=1:length(anl.Cf) plot(anl.Cf(j)+0*anl.y3_5.UyyP(j,:),anl.y3_5.UyyP(j,:)-anl.y3_5.UyyM(j,:),'r.');hold all; end
% my_legend_add(anl.lgnd);
% for j=1:length(anl.Cf) plot(anl.Cf(j)+0*anl.y7_5.UyyP(j,:),anl.y7_5.UyyP(j,:)-anl.y7_5.UyyM(j,:),'b.');hold all; end
% my_legend_add(anl.lgnd);
plot(mean(anl.y3_5.Cf(index_e,index_sg),2),mean(anl.y3_5.UyyP(index_e,index_sg),2),'ro');
%plot(mean(anl.y7_5.Cf(index_e,index_sg),2),mean(anl.y7_5.UyyP(index_e,index_sg),2),'bo');
ylabel('UyyP');





