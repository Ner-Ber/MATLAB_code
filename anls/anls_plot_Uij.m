function anls_plot_Uij(e,c)

clear c;
c1=[0.6350 0.0780 0.1840];
c2=[0 0.4470 0.7410];
for j=1:length(e.acqE.x_sg)
    if(e.acqE.x_sg(j)<100)
        c{j}=c1;
    else
        c{j}=c2;
    end
end

%Supershear
%x1=-0.01;
%x2=0.02;
x1=-50;
x2=50;

%---to cut the x axis



tOffset=00;
smt=1;
NormFactor=2;
% ----Uxy
fig1=figure;

%plot(e.acqE.tOffset-tOffset,(e.acqE.Uxy-repmat( e.pre.Uxy, length(e.acqE.tOffset(:,1)),1 ))/NormFactor,'.-');
%plot(e.acqE.tOffset-tOffset,e.acqE.Uxy-repmat( e.anl.Uxyf, length(e.acqE.tOffset(:,1)),1 ),'.-');
%plot(e.acqEc.intVdtMat-tOffset,my_smooth(e.acqEc.Uxy,smt)-repmat( e.anl.Uxyf, length(e.acqEc.intVdtMat(:,1)),1 ),'.-' );
plot(e.acqE.intVdtMat-tOffset,(my_smooth(e.acqE.Uxy,smt)-repmat( e.anl.Uxyf, length(e.acqE.intVdtMat(:,1)),1 )),'.-' );
%plot(e.acqE.tv-tOffset,my_smooth(e.acqE.Uxy,smt)-repmat( e.anl.Uxyf, length(e.acqE.intVdtMat(:,1)),1 ),'.-' );
%plot(e.acqE.intVdtMat-tOffset,my_smooth(e.acqE.Uxy,smt)-repmat( e.pre.Uxy, length(e.acqE.intVdtMat(:,1)),1 ),'.-' );
xlim([x1 x2]);
%ylim([-0.4 0.1]);
hold all;
legend([num2str(e.acqE.x_sg'), repmat('---',length(e.acqE.x_sg),1), num2str(e.acqE.y_sg')] );legend off
title([e.acqE.Date '--' e.acqE.exp '--E' num2str(e.acqE.event)]);
a=get(gca,'Children');
% if (nargin>1)
%     for j=1:length(a) set(a(end-j+1),'color',c{end-j+1}); end
% end
%-------Uxx
fig2=figure;
%plot(e.acqE.tOffset-tOffset,(e.acqE.Uxx-repmat( e.pre.Uxx, length(e.acqE.tOffset(:,1)),1 ))/NormFactor,'.-');
plot(e.acqE.intVdtMat-tOffset,my_smooth(e.acqE.Uxx,smt)-repmat( e.pre.Uxx , length(e.acqE.intVdtMat(:,1)),1 ),'.-');
xlim([x1 x2]);
%ylim([-0.9 0.05]);
hold all;
legend([num2str(e.acqE.x_sg'), repmat('---',length(e.acqE.x_sg),1), num2str(e.acqE.y_sg')] );legend off
title([e.acqE.Date '--' e.acqE.exp '--E' num2str(e.acqE.event) '- ' num2str(e.acqE.gV)]);
a=get(gca,'Children');
% if (nargin>1)
%     for j=1:length(a) set(a(end-j+1),'color',c{end-j+1}); end
% end

% % % %--For extra Uxx strain gages
% plot(e.acqE.s.tOffset_Uxx-tOffset,subtruct_norm(e.acqE.s.Uxx),'.-')
% my_legend_add(num2str(e.acqE.s.x_Uxx'));legend off;
% a=get(gca,'Children');
% if (nargin>1)
% for j=1:length(a) set(a(end-j+1),'color',c{end-j+1}); end
% end

% %-------Uyy
fig3=figure;

%plot(e.acqE.tOffset-tOffset,(e.acqE.Uyy-repmat( e.pre.Uyy, length(e.acqE.tOffset(:,1)),1 ))/NormFactor,'.-');
plot(e.acqE.intVdtMat-tOffset,my_smooth(e.acqE.Uyy,smt)-repmat( e.pre.Uyy , length(e.acqE.intVdtMat(:,1)),1 ),'.-');
xlim([x1 x2]);
%ylim([-0.1 0.3]);
hold all;
legend([num2str(e.acqE.x_sg'), repmat('---',length(e.acqE.x_sg),1), num2str(e.acqE.y_sg')] );legend off
title([e.acqE.Date '--' e.acqE.exp '--E' num2str(e.acqE.event)]);
a=get(gca,'Children');
% if (nargin>1)
%     for j=1:length(a) set(a(end-j+1),'color',c{end-j+1}); end
% end

% % %--For extra Uyy  strain gages
%  plot(e.acqE.s.tOffset_Uyy-tOffset,subtruct_norm( e.acqE.s.Uyy),'.-')
% my_legend_add(num2str(e.acqE.s.x_Uyy'));legend off
% a=get(gca,'Children');
% if (nargin>1)
% for j=1:length(a) set(a(end-j+1),'color',c{end-j+1}); end
% end

 my_createfigure([fig1 fig2 fig3]);
 close([fig1 fig2 fig3])

