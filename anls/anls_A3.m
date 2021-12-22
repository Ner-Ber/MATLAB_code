function anl=anls_A3(exper,event,slow_event)

% unlike anls_A and anls_A2 this function doesn't take average values from different
% strain gages.


%-------------
path=pwd;
path=[path(end-9:end) '\' exper '\e'];

tmp=figure;



e_Slow=[];

for j=1:length(event)
    
    tmp=find(slow_event-event(j)==0);
    if(isempty(tmp))
        kind='Rayleigh';
    else
        kind='Slow';
        e_Slow=[e_Slow event(j)];
    end
    
    %e=anls_front_in_space(exper,event(j),50,185,kind);
    e=anls_front_in_space(exper,event(j),70,160,kind);
    lgnd{j}=[path num2str(event(j))];
    
    front{j}=e.front;
    frontRaw{j}.xv=e.frontRaw.xv;
    frontRaw{j}.v=e.frontRaw.v;
    phECut{j}=e.phECut;      
    phECut{j}.lines=[];
    phECut{j}.xOffset=[];
    phECut{j}.firstline=e.phE.firstLine;
    %phECut{j}.xf_2=e.phECut.xf_2;
    %phECut{j}.xf_1=e.phECut.xf_1;
    
    
   
    
    %--------plot for check
    
    xStart=90;
    xEnd=130;
    [~,index1]=min(abs( e.phECut.frontX-xStart ));
    [~,index2]=min(abs( e.phECut.frontX-xEnd ));
    
    
    
    for k=index1:index2
        %plot(e.phECut.xOffset(k,:),(e.phECut.lines(k,:)-e.phECut.Af(k))/(1-e.phECut.Af(k)),'b.-');
        plot(e.phECut.xOffset(k,:),(e.phECut.lines(k,:)-e.front.Af(k))/(1-e.front.Af(k)),'b.-');
        %plot(e.phECut.xOffset(k,:),e.phECut.lines(k,:),'b.-');
        %figure(tmp);
%         [~,index1]=min(abs(e.phECut.x-(97-5)));
%         [~,index2]=min(abs(e.phECut.x-(97+5)));
%         index=round(linspace(index1,index2,10));
%         plot(e.phECut.xOffset(:,index),e.phECut.lines(:,index),'b.-');
         my_legend_add(num2str(e.phECut.frontX(k)'));
        hold all;
        
%         figure(21);
%         [~,index1]=min(abs(e.phECut.x-(150-5)));
%         [~,index2]=min(abs(e.phECut.x-(150+5)));
%         index=round(linspace(index1,index2,3));
%         plot(e.phECut.xOffset(:,index),e.phECut.lines(:,index),'r.-');
%         my_legend_add(num2str(e.phECut.x(index)'));
%         hold all;
        
    end
    xlim([-15 5]);
    title(num2str(event(j)));
    %pause;
    XcTmp= e.front.x1-e.front.x2;
%     XcTmp=front{j}.x1-front{j}.x2+0.25;
     Xc=mean(XcTmp(front{j}.x1>xStart&front{j}.x1<xEnd));
%         
     plot(-(Xc),0.37,'Redo');
%     plot(-front{j}.xf,0,'Redo')
%     xlim([-12.5 2.5]);
%     x=-30:0.1:30;
%     A=exp(x/Xc);
%     A(x>0)=1;
%     plot(x+0.25,A,'black.-')
%     title(num2str(event(j)));
%     hold off;
%     
     %pause on;
     pause(0.5);
    
    
    
    % plot(anl.x{j},(anl.A{j}-anl.Af{j})/(1-anl.Af{j}),'black.-');
    % xlim([-20 5]);
    %legend(num2str(e.phECut.frontX(index1:index2)) );legend off;
    
    %assignin('caller','anl',anl);
    
    
end

anl.lgnd=lgnd;
anl.event=event;
anl.event_slow=e_Slow;
anl.front=front;
anl.frontRaw=frontRaw;
anl.phECut=phECut;

