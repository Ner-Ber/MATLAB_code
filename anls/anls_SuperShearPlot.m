function anlTmp=anls_SuperShearPlot

%path1='C:\Users\IlyaS\Desktop\SuperShear\2016-05-29';
%path2={'2016-05-29-16-29-53\SubRayleigh', '2016-05-29-17-18-29\SubRayleigh' ,'2016-05-29-18-01-38\SubRayleigh', '2016-05-29-18-57-09\SubRayleigh'};%,'2016-05-29-19-31-25'
%path2={'2016-05-29-16-29-53', '2016-05-29-17-18-29' ,'2016-05-29-18-01-38', '2016-05-29-18-57-09'};%,'2016-05-29-19-31-25'

%path1='C:\Users\IlyaS\Desktop\SuperShear\2016-05-25';
%path2={'2016-05-25-18-49-33\SubRayleigh', '2016-05-25-19-39-20\SubRayleigh'};
%path2={'2016-05-25-18-49-33', '2016-05-25-19-39-20'};%,'2016-05-29-19-31-25'

%path1='C:\Users\IlyaS\Desktop\SuperShear\2016-05-27';
%path2={'2016-05-27-12-12-30\SubRayleigh', '2016-05-27-12-49-25\SubRayleigh'};%,'2016-05-29-19-31-25'
%path2={'2016-05-27-12-12-30', '2016-05-27-12-49-25'};

path1='c:\frics\2013-01-29';
path2={'18-33-32\SubRayleigh'};


% x3_5Index=[10,12,14,15];
% x7_5Index=[11,13];

x3_5Index=[9:12];


% for j=1:length(path2)
%     
%     load([path1 '\' path2{j} '\anl.mat']) ;
%     anl.y3_5.Uxy0=anl.pre.Uxy(:,x3_5Index);
%     anl.y3_5.Uxyf=anl.y3_5.Uxy0-anl.y3_5.Uxy0_f;
%     
% %    anl.y7_5.Uxy0=anl.pre.Uxy(:,x7_5Index);
%  %  anl.y7_5.Uxyf=anl.y7_5.Uxy0-anl.y7_5.Uxy0_f;
%     
%     anlTmp(j)=anl;
%     
% end

anl=anlTmp;


%------------------Plot y=3.5mm Uxy_r Vs Uxy_0
% figN=figure;
% %----Compare different experiments
% for j=1:length(anl)
%     for k=1:length(anl(j).y3_5.x_sg)
%         figure(figN+k-1);
%         plot(anl(j).y3_5.Uxy0(:,k),anl(j).y3_5.Uxyf(:,k),'o');
%         a=fit(anl(j).y3_5.Uxy0(:,k),anl(j).y3_5.Uxyf(:,k),'poly1');
%         a=confint(a);
%         p1=mean(a(:,1));
%         dp1=mean( a(:,1) ) - a(1,1);
%         hold all;
%         xlabel('Uxy_0');
%         ylabel('Uxy_r');
%         my_legend_add([anl(j).lgnd{1} '--p1=' num2str(p1) '+-' num2str(dp1) ]);
%         title(['y=3.5,x=' num2str(anl(j).y3_5.x_sg(k))]);
%     end
% end

% for k=1:length(anl(j).y3_5.x_sg)
%     figure(figN+k-1);
%     fig=my_get_axis;
%     x=fig.x{1};
%     y=fig.y{1};
%     for l=2:length(fig.x)
%         x=[x fig.x{l}];
%         y=[y fig.y{l}];
%     end
%     plot(x,y,'blacko');
%     a=fit(x',y','poly1');
%     plot(a);
%     a=confint(a);
%     title([fig.title ',p1=' num2str(a(:,1)') ',p2=' num2str(a(:,2)') ]);
%     xlim([1.2 1.7]);
%     ylim([1.1 1.5]);
% end
%
% %----------Plot y=7.5mm
% for j=1:length(anl)
%
%     for k=1:length(anl(j).y7_5.x_sg)
%         figure(29+k);
%         plot(anl(j).y7_5.Uxy0(:,k),anl(j).y7_5.Uxyf(:,k),'o');
%         hold all;
%         my_legend_add(anl(j).lgnd(1));
%         title(['y=7.5,x=' num2str(anl(j).y7_5.x_sg(k))]);
%     end
%
% end
%
% for k=1:length(anl(j).y7_5.x_sg)
%     figure(29+k);
%     fig=my_get_axis;
%     x=fig.x{1};
%     y=fig.y{1};
%     for l=2:length(fig.x)
%         x=[x fig.x{l}];
%         y=[y fig.y{l}];
%     end
%     plot(x,y,'blacko');
%     a=fit(x',y','poly1');
%     plot(a);
%     a=confint(a);
%     title([fig.title ',p1=' num2str(a(:,1)') ',p2=' num2str(a(:,2)') ]);
%     xlim([1.2 1.7]);
%     ylim([1.1 1.5]);
% end


% %----Compare different experiments and plot tau_r Vs Uxx*Cf figure for each sg location
mu_0=2.1;
mu_inf=1.2;
figN=figure;
for j=1:length(anl)
    for k=1:length(anl(j).y3_5.x_sg)
        subplot(2,2,k);
        tau_r=2*mu_0*(anl(j).y3_5.Uxyf(:,k)-anl(j).y3_5.Uxy0(:,k)*(1-mu_inf/mu_0));
        v=anl(j).y3_5.UxxP(:,k)*1e-3.*anl(j).Cf';
        plot(v,tau_r,'o');
        xlabel('Uxx*Cf(m/s)');
        ylabel('tau_r(Mpa)');
        hold all;
        my_legend_add(anl(j).lgnd(1));
        title(['y=3.5,x=' num2str(anl(j).y3_5.x_sg(k)) '\mu_inf=' num2str(mu_inf) ',mu_0= ' num2str(mu_0) ]);
    end
end




%-------Join different experiments and plot --- y=3.5mm
% Uxy0=anl(1).y3_5.Uxy0;
% Uxyf=anl(1).y3_5.Uxyf;
% UxxP=anl(1).y3_5.UxxP;
% Cf=anl(1).Cf;
%
% for j=2:length(anl)
%     Uxy0=[Uxy0 ; anl(j).y3_5.Uxy0];
%     Uxyf=[Uxyf ;anl(j).y3_5.Uxyf];
%     UxxP=[UxxP ;anl(j).y3_5.UxxP];
%     Cf=[Cf anl(j).Cf];
% end
%
% UxxPMean=mean(UxxP,2);
%
% figure(figN);
% %plot(Uxy0,Uxyf,'o');
% tau_r=5.65*(Uxyf-Uxy0*(1-3/5.65));
% v=UxxPMean*1e-3.*Cf';
% plot(v,tau_r,'o');
%
% hold all;

%------- plot individual experiments --- y=3.5mm

% for j=1:length(anl)
% figure(figN+j);
%     Uxy0=anl(j).y3_5.Uxy0;
%     Uxyf=anl(j).y3_5.Uxyf;
%     UxxP=anl(j).y3_5.UxxP;
%     Cf=anl(j).Cf';
%     v=UxxP*1e-3.*repmat(Cf,1,length(UxxP(1,:)));
% tau_r=5.65*(Uxyf-Uxy0*(1-3/5.65));
% plot(v,tau_r,'o');
% title(anl(j).lgnd(1));
% end
%
% hold all;
%%---------------


%-------Join different experiments fit and plot --- 7=3.5mm
% Uxy0=anl(1).y7_5.Uxy0;
% Uxyf=anl(1).y7_5.Uxyf;
% UxxP=anl(1).y7_5.UxxP;
% Cf=anl(1).Cf;
%
% for j=2:length(anl)
%     Uxy0=[Uxy0 ; anl(j).y7_5.Uxy0];
%     Uxyf=[Uxyf ;anl(j).y7_5.Uxyf];
%     UxxP=[UxxP ;anl(j).y7_5.UxxP];
%     Cf=[Cf anl(j).Cf];
% end
%
% figure(figN);
% tau_r=5.65*(Uxyf-Uxy0*(1-3/5.65));
% v=UxxP*1e-3.*repmat(Cf',1,length(UxxP(1,:)));
% plot(v,tau_r,'o');
%
% %plot(Uxy0,Uxyf,'o');
%
% hold all;

% for j=1:length(Uxy0(1,:))
% a=fit(Uxy0(:,j),Uxyf(:,j),'poly1');
% plot(a);
% a=confint(a);
% my_legend_add( [',p1=' num2str(a(:,1)') ',p2=' num2str(a(:,2)') ] );
% end
%
% xlim([1.2 1.7]);
% ylim([1.1 1.5]);






