function anl=anls_A(exper,event)


% x3_5Index=[10,12,14,15];
% x7_5Index=[11,13];
% xStart=110;
% xEnd=130;

x3_5Index=[9:12];
x7_5Index=[ ];
xStart=105;
xEnd=130;


%-------------
path=pwd;
path=[path(end-9:end) '\' exper '\e'];

tmp=figure;


for j=1:length(event)
    lgnd=[path num2str(event(j))];
    e=anls_front_in_space(exper,event(j),30,185);
    
    anl.lgnd{j}=lgnd;
    anl.event(j)=event(j);
    
    %------calc mean Cf
    [~,index1]=min(abs(e.frontRaw.xv-xStart ));
    [~,index2]=min(abs(e.frontRaw.xv-xEnd ));
    anl.Cf(j)=mean(e.frontRaw.v(index1:index2)) ;
    
    %---
    [~,index1]=min(abs( e.phECut.frontX-xStart ));
    [~,index2]=min(abs( e.phECut.frontX-xEnd ));
    
    [anl.x{j},anl.A{j}]=averageVec(e.phECut.xOffset(index1:index2,:)',e.phECut.lines(index1:index2,:)');
    anl.Af(j)=mean(e.phECut.Af(index1:index2));
    anl.AfStd(j)=std(e.phECut.Af(index1:index2));
    anl.xf(j)=e.front.xf;
    
    
    Xc=e.front.x1-e.front.x2+0.25;
    anl.Xc(j)=mean(Xc(e.front.x1>xStart&e.front.x1<xEnd));
    anl.XcStd(j)=std(Xc(e.front.x1>xStart&e.front.x1<xEnd));
    
    anl.UxxP3_5(j)=mean(e.anl.UxxP(x3_5Index));
    anl.UxxP3_5Std(j)=std(e.anl.UxxP(x3_5Index));
    
    %anl.UxxP7_5(j)=mean(e.anl.UxxP(x7_5Index));
    %anl.UxxP7_5Std(j)=std(e.anl.UxxP(x7_5Index));
    
    anl.UyyP3_5(j)=mean(e.anl.UyyP(x3_5Index));
    anl.UyyP3_5Std(j)=std(e.anl.UyyP(x3_5Index));
    
    %anl.UyyP7_5(j)=mean(e.anl.UyyP(x7_5Index));
    %anl.UyyP7_5Std(j)=std(e.anl.UyyP(x7_5Index));
    
    %--------plot for check
    subplot(1,3,1);
    for k=index1:index2
        plot(e.phECut.xOffset(k,:),(e.phECut.lines(k,:)-e.phECut.Af(k))/(1-e.phECut.Af(k)),'b.-');
        hold all;
    end
    plot(-anl.Xc(j),0.37,'Redo');
    plot(-anl.xf(j),0,'Redo')
    xlim([-12.5 2.5]);
    x=-30:0.1:30;
    A=exp(x/anl.Xc(j));
    A(x>0)=1;
    plot(x+0.25,A,'black.-')
    title(num2str(event(j)));
    hold off;
    
    subplot(1,3,2);
    plot(e.acqE.t,e.acqE.Uxx(:,[x3_5Index x7_5Index])-repmat(e.pre.Uxx([x3_5Index x7_5Index]),length(e.acqE.t),1),'.-');
    hold all;
    plot(e.anl.UxxPt([x3_5Index x7_5Index]),-e.anl.UxxP([x3_5Index x7_5Index]),'blacko');
    xlim([e.anl.UxxPt(x3_5Index(1))-0.03 e.anl.UxxPt(x3_5Index(1))+0.03]);
    hold off
    
    subplot(1,3,3);
    plot(e.acqE.t,e.acqE.Uyy(:,[x3_5Index x7_5Index])-repmat(e.pre.Uyy([x3_5Index x7_5Index]),length(e.acqE.t),1),'.-');
    hold all;
    plot(e.anl.UyyPt([x3_5Index x7_5Index]),e.anl.UyyP([x3_5Index x7_5Index]),'blacko');
    xlim([e.anl.UxxPt(x3_5Index(1))-0.03 e.anl.UxxPt(x3_5Index(1))+0.03]);
    hold off
    
    pause on;
    pause(0.5);

    
    
    % plot(anl.x{j},(anl.A{j}-anl.Af{j})/(1-anl.Af{j}),'black.-');
    % xlim([-20 5]);
    %legend(num2str(e.phECut.frontX(index1:index2)) );legend off;
    
    %assignin('caller','anl',anl);
    
    
    
    
end