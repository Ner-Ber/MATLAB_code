function anlEqm=anls_SuperShear_Eqm(exper,event)

%Very similar to anls_SuperShear but analyzes only Cf Vs. x and Uxy_0-Uxy_r

path=pwd;
path=[path(end-9:end) '\' exper '\e'];

tmp=figure;

for j=1:length(event)
    j/length(event)
    lgnd=[path num2str(event(j))];
    e=anls_front_in_space(exper,event(j),40,180,'Super');
    
    x3_5Index=1:length(e.acqE.x_sg);
    x7_5Index=1:length(e.acqE.x_sg);
    x3_5Index=x3_5Index(e.acqE.y_sg(x3_5Index)<5);
    x7_5Index=x7_5Index(e.acqE.y_sg(x7_5Index)>5);
    
    anl.lgnd{j}=lgnd;
    anl.e(j)=event(j);
    anl.pre.x_sg=e.acqE.x_sg;
    anl.pre.Uxy(j,:)=e.pre.Uxy;
    anl.pre.Uyy(j,:)=e.pre.Uyy;
    
    anl.frontRaw{j}.xv=e.frontRaw.xv;
    anl.frontRaw{j}.v=e.frontRaw.v;
    
    anl.y3_5.x_sg=e.acqE.x_sg(x3_5Index);
    anl.y3_5.Uxy0_f(j,:)=e.pre.Uxy(x3_5Index)-e.anl.Uxyf(x3_5Index);
    
    anl.y7_5.x_sg=e.acqE.x_sg(x7_5Index);
    anl.y7_5.Uxy0_f(j,:)=e.pre.Uxy(x7_5Index)-e.anl.Uxyf(x7_5Index);
    
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
        anl.y3_5.Uxy0_f(j,:)=e.pre.Uxy(x3_5Index) - e.anl.Uxyf(x3_5Index) - y(1:length(x3_5Index));
        anl.y7_5.Uxy0_f(j,:)=e.pre.Uxy(x7_5Index) - e.anl.Uxyf(x7_5Index) - y(length(x3_5Index)+1:end);
    end
    
     
end

anlEqm=anl;

