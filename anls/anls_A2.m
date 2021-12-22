function anl=anls_A2(exper,event)

% unlike anls_A this function doesn't take average values from different
% strain gages.


%-------------
path=pwd;
path=[path(end-9:end) '\' exper '\e'];

tmp=figure;

Cf=zeros(length(event),length(x3_5Index));
A0=Cf;
Af=Cf;
AfStd=Cf;
xf=Cf;
Xc=Cf;
XcStd=Cf;


for j=1:length(event)
    
    e=anls_front_in_space(exper,event(j),40,170,'Super');
    lgnd{j}=[path num2str(event(j))];


    x3_5Index=5:length(e.acqE.x_sg);
    x7_5Index=1:length(e.acqE.x_sg);
    x3_5Index=x3_5Index(e.acqE.y_sg(x3_5Index)<5);
    x7_5Index=x7_5Index(e.acqE.y_sg(x7_5Index)>5);
    
    %----- Cf
    
    for k=1:length(x3_5Index)
        
        xStart=e.acqE.x_sg(x3_5Index(k))-5;
        xEnd=e.acqE.x_sg(x3_5Index(k))+5;
        Cf(j,k)=e.frontSg.v(x3_5Index(k));
        %---For each strain gage take average over +-5mm
        A0(j,k)=mean(e.phE.firstLine((e.phE.x>xStart&e.phE.x<xEnd)));
        Af(j,k)=mean(e.phECut.Af((e.front.x1>xStart&e.front.x1<xEnd)));
        AfStd(j,k)=std(e.phECut.Af((e.front.x1>xStart&e.front.x1<xEnd)));
        xf(j)=e.front.xf;
        
        XcTmp=e.front.x1-e.front.x2+0.25;
        Xc(j,k)=mean(XcTmp(e.front.x1>xStart&e.front.x1<xEnd));
        XcStd(j,k)=std(XcTmp(e.front.x1>xStart&e.front.x1<xEnd));
        
    end
    %--------plot for check
    k_chk=3;
    xStart=e.acqE.x_sg(x3_5Index(k_chk))-5;
    xEnd=e.acqE.x_sg(x3_5Index(k_chk))+5;
    [~,index1]=min(abs( e.phECut.frontX-xStart ));
    [~,index2]=min(abs( e.phECut.frontX-xEnd ));
    
    for k=index1:index2
        plot(e.phECut.xOffset(k,:),(e.phECut.lines(k,:)-e.phECut.Af(k))/(1-e.phECut.Af(k)),'b.-');
        hold all;
    end
    plot(-(Xc(j,k_chk)-0.25),0.37,'Redo');
    plot(-xf(j),0,'Redo')
    xlim([-12.5 2.5]);
    x=-30:0.1:30;
    A=exp(x/Xc(j,k_chk));
    A(x>0)=1;
    plot(x+0.25,A,'black.-')
    title(num2str(event(j)));
    hold off;
    
    pause on;
    pause(0.5);
    
    
    
    % plot(anl.x{j},(anl.A{j}-anl.Af{j})/(1-anl.Af{j}),'black.-');
    % xlim([-20 5]);
    %legend(num2str(e.phECut.frontX(index1:index2)) );legend off;
    
    %assignin('caller','anl',anl);
    
    
end

anl.lgnd=lgnd;
anl.event=event;
anl.Cf=Cf;
anl.A0=A0;
anl.Af=Af;
anl.AfStd=AfStd;
anl.xf=xf;
anl.Xc=Xc;
anl.XcStd=XcStd;
anl.x_sg=e.acqE.x_sg(x3_5Index);