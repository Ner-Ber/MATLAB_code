function anl=anls_EquationOfMotion(exper,event,slow_event)
%The functions helps to get Uxy0-Uxy_r, UxxP and UyyP for equation of
%mation analysys. Use the function anls_EquationOfMotion_Plot(anl) to plot.

check=0;
path=pwd;
path=[path(end-9:end) '\' exper '\e'];

tmp_fig=figure;

for j=1:length(event)
    j/length(event)
    lgnd=[path num2str(event(j))];
    
    tmp=find(slow_event-event(j)==0);
    if(isempty(tmp))
        kind='Rayleigh';
    else
        kind='Slow';
    end
    
    e=anls_front_in_space(exper,event(j),50,180,kind);%45
    
    x3_5Index=1:length(e.acqE.x_sg);
    %     x7_5Index=1:length(e.acqE.x_sg);
    x3_5Index=x3_5Index(e.acqE.y_sg(x3_5Index)<5);
    
    %x3_5Index=x3_5Index(6:end-5);%----------Not all sg
    %x3_5Index=x3_5Index(1:end);%----------
    
    %     x7_5Index=x7_5Index(e.acqE.y_sg(x7_5Index)>5);
    %
    %          x3_5Index=1:length(e.acqE.x_sg);
    x7_5Index=[];
    
    anl.lgnd{j}=lgnd;
    anl.e(j)=event(j);
    
    anl.frontRaw{j}.xv=e.frontRaw.xv;
    anl.frontRaw{j}.v=e.frontRaw.v;
    
    anl.y3_5.x_sg=e.acqE.x_sg(x3_5Index);
    anl.y3_5.Uxy0_f(j,:)=e.pre.Uxy(x3_5Index)-e.anl.Uxyf(x3_5Index);
    anl.y3_5.Uxy0(j,:)=e.pre.Uxy(x3_5Index);
    anl.y3_5.Uxx0(j,:)=e.pre.Uxx(x3_5Index);
    anl.y3_5.Uyy0(j,:)=e.pre.Uyy(x3_5Index);
    anl.y3_5.UxxP(j,:)=e.anl.UxxP(x3_5Index);
    anl.y3_5.UyyP(j,:)=e.anl.UyyP(x3_5Index);
    anl.y3_5.Uxyf(j,:)=e.anl.Uxyf(x3_5Index);
    anl.y3_5.Uxxf(j,:)=e.anl.Uxxf(x3_5Index);
    anl.y3_5.xf=e.anl.xf;
    
    %---extra locations
    anl.y3_5.xf2=e.anl.xf2;
    anl.y3_5.Uxxf2(j,:)=e.anl.Uxxf2(x3_5Index);
    anl.y3_5.xf3=e.anl.xf3;
    anl.y3_5.Uxxf3(j,:)=e.anl.Uxxf3(x3_5Index);
    anl.y3_5.Uxyf3(j,:)=e.anl.Uxyf3(x3_5Index);
     anl.y3_5.Uxy0_f3(j,:)=e.pre.Uxy(x3_5Index)-e.anl.Uxyf3(x3_5Index);    
    %     anl.y7_5.x_sg=e.acqE.x_sg(x7_5Index);
    %     anl.y7_5.Uxy0_f(j,:)=e.pre.Uxy(x7_5Index)-e.anl.Uxyf(x7_5Index);
    %     anl.y7_5.Uxy0(j,:)=e.pre.Uxy(x7_5Index);
    %     anl.y7_5.Uxyf(j,:)=e.anl.Uxyf(x7_5Index);
    %     anl.y7_5.UxxP(j,:)=e.anl.UxxP(x7_5Index);
    %     anl.y7_5.Uxy0(j,:)=e.pre.Uxy(x7_5Index);
    %     anl.y7_5.UyyP(j,:)=e.anl.UyyP(x7_5Index);
    
    if(check ==1)
        %------------check and fix
        figure(tmp_fig);
        
        plot(e.acqE.intVdtMat(:,[x3_5Index x7_5Index]),e.acqE.Uxy(:,[x3_5Index x7_5Index])-repmat( e.anl.Uxyf3([x3_5Index x7_5Index]), length(e.acqE.intVdtMat(:,1)),1 ),'.-' );
        xlim([-50 50]);
        pause(0.2);
        
        title(num2str(anl.e(j)))
        state=input('To fix press y :','s');
        %state='n';
        if strcmp(state,'y')
            [x,y]=my_ginput;
            anl.y3_5.Uxy0_f(j,:)=e.pre.Uxy(x3_5Index) - e.anl.Uxyf3([x3_5Index ]) - y(1:length(x3_5Index));
            %      anl.y7_5.Uxy0_f(j,:)=e.pre.Uxy(x7_5Index) - e.anl.Uxyf([x7_5Index ]) - y(length(x3_5Index)+1:end);
        end
        
        %%-------measure Uxxf
        figure(tmp_fig);
        plot(e.acqE.intVdtMat(:,[x3_5Index x7_5Index]),e.acqE.Uxx(:,[x3_5Index x7_5Index])-repmat( e.pre.Uxx([x3_5Index x7_5Index]), length(e.acqE.intVdtMat(:,1)),1 ),'.-' );
        xlim([-50 50]);
        pause(0.2);
        
        hold all;
        plot(e.anl.UxxP([x3_5Index x7_5Index])*0+anl.y3_5.xf,-e.anl.Uxxf([x3_5Index x7_5Index]),'blacko');
        hold off;
        title(num2str(anl.e(j)))
        pause(0.5)
%         state=input('To fix press y :','s');
%         if strcmp(state,'y')
%             [x,y]=my_ginput;
%             anl.y3_5.Uxxf(j,:)=-y(1:length(x3_5Index));
%             anl.y7_5.Uxxf(j,:)=-y(length(x3_5Index)+1:end);
%         end

%         %-------measure UxxP
%         figure(tmp_fig);
%         plot(e.acqE.intVdtMat(:,[x3_5Index x7_5Index]),e.acqE.Uxx(:,[x3_5Index x7_5Index])-repmat( e.pre.Uxx([x3_5Index x7_5Index]), length(e.acqE.intVdtMat(:,1)),1 ),'.-' );
%         xlim([-50 50]);
%         pause(0.2);
%         
%         hold all;
%         plot(e.anl.UxxPx([x3_5Index x7_5Index]),-e.anl.UxxP([x3_5Index x7_5Index]),'blacko');
%         hold off;
%         title(num2str(anl.e(j)))
%         
%         state=input('To fix press y :','s');
%         if strcmp(state,'y')
%             [x,y]=my_ginput;
%             anl.y3_5.UxxP(j,:)=-y(1:length(x3_5Index));
%             anl.y7_5.UxxP(j,:)=-y(length(x3_5Index)+1:end);
%         end
        
%         %-------measure UyyP
%         figure(tmp_fig);
%         plot(e.acqE.intVdtMat(:,[x3_5Index x7_5Index]),e.acqE.Uyy(:,[x3_5Index x7_5Index])-repmat( e.pre.Uyy([x3_5Index x7_5Index]), length(e.acqE.intVdtMat(:,1)),1 ),'.-' );
%         xlim([-50 50]);
%         pause(0.2);
%         
%         hold all;
%         plot(e.anl.UyyPx([x3_5Index x7_5Index]),e.anl.UyyP([x3_5Index x7_5Index]),'blacko');
%         hold off;
%         title(num2str(anl.e(j)))
%         
%         state=input('To fix press y :','s');
%         if strcmp(state,'y')
%             [x,y]=my_ginput;
%             anl.y3_5.UyyP(j,:)=y(1:length(x3_5Index));
%             anl.y7_5.UyyP(j,:)=y(length(x3_5Index)+1:end);
%         end
    end
end


