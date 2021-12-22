function anl=anls_SuperShear(exper,event)

%x3_5Index=[11,13,15,16];
%x7_5Index=[12,14];

%x3_5Index=[8,12];
%x7_5Index=(9:11);

%x3_5Index=[10,12,14,15];
%x7_5Index=[11,13];

%x3_5Index=[7:19];
%x7_5Index=[];

path=pwd;
path=[path(end-9:end) '\' exper '\e'];

tmp=figure;

for j=1:length(event)
    j/length(event)
    lgnd=[path num2str(event(j))];
    e=anls_front_in_space(exper,event(j),50,170,'Super');
    
    x3_5Index=5:length(e.acqE.x_sg);
    x7_5Index=1:length(e.acqE.x_sg);
    x3_5Index=x3_5Index(e.acqE.y_sg(x3_5Index)<5);
    x7_5Index=x7_5Index(e.acqE.y_sg(x7_5Index)>5);
    
    anl.lgnd{j}=lgnd;
    anl.e(j)=event(j);
    anl.pre.x_sg=e.acqE.x_sg;
    anl.pre.Uxy(j,:)=e.pre.Uxy;
    anl.pre.Uyy(j,:)=e.pre.Uyy;
    
    %anl.Cf(j)=mean(e.frontSg.v(min([x3_5Index x7_5Index]):max([x3_5Index x7_5Index])));
%     [~,index1]=min(abs(e.frontRaw.xv-e.acqE.x_sg(x3_5Index(1)) ));
%     [~,index2]=min(abs(e.frontRaw.xv-e.acqE.x_sg(x3_5Index(end))));
%     anl.Cf(j)=mean(e.frontRaw.v(index1:index2)) ;
    
    anl.frontRaw{j}.xv=e.frontRaw.xv;
    anl.frontRaw{j}.v=e.frontRaw.v;
    
    anl.y3_5.x_sg=e.acqE.x_sg(x3_5Index);
    anl.y3_5.Uxy0_f(j,:)=e.pre.Uxy(x3_5Index)-e.anl.Uxyf(x3_5Index);
    anl.y3_5.UxxP(j,:)=e.anl.UxxP(x3_5Index);
    anl.y3_5.UyyP(j,:)=e.anl.UyyP(x3_5Index);
    anl.y3_5.Uyy1(j,:)=e.anl.Uyy1(x3_5Index);    
    
     anl.y7_5.x_sg=e.acqE.x_sg(x7_5Index);
     anl.y7_5.Uxy0_f(j,:)=e.pre.Uxy(x7_5Index)-e.anl.Uxyf(x7_5Index);
     anl.y7_5.UxxP(j,:)=e.anl.UxxP(x7_5Index);
     anl.y7_5.UyyP(j,:)=e.anl.UyyP(x7_5Index);
     anl.y7_5.Uyy1(j,:)=e.anl.Uyy1(x7_5Index);
    
%------------check and fix
   figure(tmp);
    
   %plot(e.acqEc.intVdtMat(:,[x3_5Index x7_5Index]),e.acqEc.Uxy(:,[x3_5Index x7_5Index])-repmat( e.pre.Uxy([x3_5Index x7_5Index]), length(e.acqEc.intVdtMat(:,1)),1 ),'.-' );
   plot(e.acqE.intVdtMat(:,[x3_5Index x7_5Index]),e.acqE.Uxy(:,[x3_5Index x7_5Index])-repmat( e.anl.Uxyf([x3_5Index x7_5Index]), length(e.acqE.intVdtMat(:,1)),1 ),'.-' );
   xlim([-70 50]);
    
   % plot(e.acqE.t,e.acqE.Uxy(:,[x3_5Index x7_5Index])-repmat(e.anl.Uxyf([x3_5Index x7_5Index]),length(e.acqE.t),1),'.-');
    %xlim([e.anl.UxxPt(x3_5Index(1))-0.03 e.anl.UxxPt(x3_5Index(1))+0.1]);
    
    title(num2str(anl.e(j)))
    pause(0.1);

    state=input('To fix press y :','s');
    
    if strcmp(state,'y')
        [x,y]=my_ginput;
        anl.y3_5.Uxy0_f(j,:)=e.pre.Uxy(x3_5Index) - e.anl.Uxyf([x3_5Index ]) - y(1:length(x3_5Index));
        anl.y7_5.Uxy0_f(j,:)=e.pre.Uxy(x7_5Index) - e.anl.Uxyf([x7_5Index ]) - y(length(x3_5Index)+1:end);
    end
    
  %---------- measure UxxP
    figure(tmp);
    plot(e.acqE.t,e.acqE.Uxx(:,[x3_5Index x7_5Index])-repmat(e.pre.Uxx([x3_5Index x7_5Index]),length(e.acqE.t),1),'.-');
    dt=20/anl.Cf(j);
    xlim([e.anl.UxxPt(x3_5Index(1))-dt e.anl.UxxPt(x3_5Index(end))+dt]);
       
    hold all;
    plot(e.anl.UxxPt([x3_5Index x7_5Index]),-e.anl.UxxP([x3_5Index x7_5Index]),'blacko');
    title([num2str(anl.e(j)) ' - Cf=' num2str(anl.Cf(j)) ]);
    hold off
    pause(0.1);
    
    state=input('To fix press y :','s');
    if strcmp(state,'y')
        [x,y]=my_ginput;
        anl.y3_5.UxxP(j,:)=-y(1:length(x3_5Index));
        anl.y7_5.UxxP(j,:)=-y(length(x3_5Index)+1:end-1);
    end
    
    %--------- Measure Mach cone in Uxx
    figure(tmp);
    plot(e.acqE.t,e.acqE.Uxx(:,[x3_5Index x7_5Index])-repmat(e.pre.Uxx([x3_5Index x7_5Index]),length(e.acqE.t),1),'.-');
    dt=20/anl.Cf(j);
    xlim([e.anl.UxxPt(x3_5Index(1))-dt e.anl.UxxPt(x3_5Index(end))+dt]);
    title([num2str(anl.e(j)) ' - Cf=' num2str(anl.Cf(j)) ]);    hold off
    display('Find the mach cone ');
    %if strcmp(state,'y')
        [x,y]=my_ginput;
        anl.y3_5.Uxx1(j,:)=-y(1:length(x3_5Index));
        anl.y7_5.Uxx1(j,:)=-y(length(x3_5Index)+1:end);
    %end
    
 % ---------------------- measure UyyP
    figure(tmp);
    plot(e.acqE.t,e.acqE.Uyy(:,[x3_5Index x7_5Index])-repmat(e.pre.Uyy([x3_5Index x7_5Index]),length(e.acqE.t),1),'.-');
    hold all;
    plot(e.anl.UyyPt([x3_5Index x7_5Index]),e.anl.UyyP([x3_5Index x7_5Index]),'blacko');
    xlim([e.anl.UxxPt(x3_5Index(1))-dt e.anl.UxxPt(x3_5Index(end))+dt]);
    title([num2str(anl.e(j)) ' - Cf=' num2str(anl.Cf(j)) ]);
    hold off;
    pause(0.1);

    state=input('To fix press y :','s');
    if strcmp(state,'y')
        [x,y]=my_ginput;
        anl.y3_5.UyyP(j,:)=y(1:length(x3_5Index));
        anl.y7_5.UyyP(j,:)=y(length(x3_5Index)+1:end-1);
    end
    
    %--------fon mach cone in Uyy
    figure(tmp);
    plot(e.acqE.t,e.acqE.Uyy(:,[x3_5Index x7_5Index])-repmat(e.pre.Uyy([x3_5Index x7_5Index]),length(e.acqE.t),1),'.-');
    hold all;
    plot(e.anl.Uyy1t([x3_5Index x7_5Index]),e.anl.Uyy1([x3_5Index x7_5Index]),'blacko');
    xlim([e.anl.UxxPt(x3_5Index(1))-0.03 e.anl.UxxPt(x3_5Index(end))+0.03]);
    title([num2str(anl.e(j)) ' - Cf=' num2str(anl.Cf(j)) ]);
    hold off;
    pause(0.1);

    state=input('To fix  mach cone press y:','s');
    if strcmp(state,'y')
        [x,y]=my_ginput;
        anl.y3_5.Uyy1(j,:)=y(1:length(x3_5Index));
        anl.y7_5.Uyy1(j,:)=y(length(x3_5Index)+1:end-1);
    end
    
    %field_names=fieldnames(anl);
    %for j=1:length(field_names)
        %assignin('caller','anl',anl); 
     %end
    
 end


