function [e]=mixVK_anls_front_in_space(exper,event,tStart,tEnd,Xstart,Xend)

acqE=mixVK_acq132_event_get_data(exper,event,'start','end',1,'Uxx','Uxy','Uyy','Sxx','Sxy','Syy','x_sg');
phE=phantomGetLines(exper,event,tStart-2,'min',tEnd+2,1000,1,'',2:7);
phE.firstLine=mean(phE.lines(1:20,:),1);
A=get_A_at_x(phE,acqE(1).x_sg,1);

%----detecte and smooth the front

[frontRaw front]=mixVK_calc_frontXT_from_XAN(phE,tStart,tEnd,3,Xstart,Xend);

%-calc front v (set the smooth parameters at the function

[frontRaw.xv frontRaw.v]=calc_front_v(frontRaw.x,frontRaw.t);

[front.xv front.v]=calc_front_v(front.xRaw,front.t);


%----------
x=frontRaw.x; %another option is by front.x1
t=frontRaw.t;

%  x=front.x; %another option is by front.x1
%  t=front.t;

%----------cutted data (in time) of phE
%phECut contains only relevent time window and also contans 'xOffset' where
%x_tip is subtructed.
%To plot use transpose first : plot(phECut.xOffset',phECut.lines','.-')

for j=1:length(t)
  [~, index_vec(j)]=min(abs(phE.t-t(j)));
end
%  [~,jS]=min(abs(phE.t-t(1)));
%  [~,jE]=min(abs(phE.t-t(end)));

phECut.t=phE.t(index_vec);
phECut.lines=phE.lines(index_vec,:)./repmat(phE.firstLine,length(phECut.t),1);
phECut.frontX=x';
phECut.x=phE.x;
phECut.xOffset=repmat(phE.x,length(phECut.t),1)-repmat(phECut.frontX,1,length(phECut.x));

%----front at sg locations
SgNumber=length(acqE.x_sg);

for j=1:SgNumber [~,index(j)]=min(abs(frontRaw.xv-acqE.x_sg(j))); end
frontSg.v=frontRaw.v(index);

for j=1:SgNumber [~,index(j)]=min(abs(front.xRaw-acqE.x_sg(j))); end
frontSg.Af=front.Af(index);

% %---simple way to calc front.t at x_sg
% SgNumber=length(acqE.x_sg);
% for j=1:SgNumber [~,index(j)]=min(abs(x-acqE.x_sg(j))); end
% frontSg.t=t(index);
% frontSg.x=x(index);

%---Better way to calc front.t at x_sg
for j=1:SgNumber
    [~,index(j)]=min(abs(x-acqE.x_sg(j)));
end
frontSg.t=t(index)+(acqE.x_sg-x(index))./frontSg.v;
frontSg.x=acqE.x_sg;

%---x=V*t
for j=1:SgNumber acqE.tv(:,j)=-(acqE.t-frontSg.t(j))*frontSg.v(j); end
for j=1:SgNumber A.tv(:,j)=-(A.t-frontSg.t(j))*frontSg.v(j); end

acqE.tOffset=repmat(acqE.t,1,length(acqE.x_sg))-repmat(frontSg.t,length(acqE.t),1);

%--- Do better then constant rupture front velocity (Integrate Vdt)

acqE.intVdt=frontX_interp(x,t,0.9999999,acqE.t);%0.999999
[~,acqE.intVdt_i1]=min(abs(acqE.t-tStart));
[~,acqE.intVdt_i2]=min(abs(acqE.t-tEnd));
for j=1:SgNumber acqE.intVdtMat(:,j)=-acqE.intVdt+frontSg.x(j); end

%acqE=structCutTime(acqE,tStart,tEnd);

A.intVdt=frontX_interp(x,t,0.9999999,A.t);
[~,A.intVdt_i1]=min(abs(A.t-tStart));
[~,A.intVdt_i2]=min(abs(A.t-tEnd));

for j=1:SgNumber A.intVdtMat(:,j)=-A.intVdt+frontSg.x(j); end
%A=structCutTime(A,tStart,tEnd);
%--- Cut the data for IntVdt
acqEc=structCutTime(acqE,acqE.t(acqE.intVdt_i1),acqE.t(acqE.intVdt_i2));
Ac=structCutTime(A,A.t(A.intVdt_i1),A.t(A.intVdt_i2));

%find Uxyf
for j=1:length(acqE.x_sg)
    [~,index1]=min(abs(acqEc.intVdtMat(:,j)--35));
    [~,index2]=min(abs(acqEc.intVdtMat(:,j)--45));
    e.anl.Uxyf(j)=mean(acqEc.Uxy(index1:index2,j),1);
end

e.frontRaw=frontRaw;
e.frontSg=frontSg;
e.front=front;
e.acqE=acqE;
e.acqEc=acqEc;
e.phE=phE;
e.phECut=phECut;
e.A=A;
e.Ac=Ac;
e.pre=mixVK_anls_pre_conditions(exper,event,-3,2:7,1);
e.anl.tStart=tStart;
e.anl.tEnd=tEnd;
e.anl.Xstart=Xstart;
e.anl.Xend=Xend;

V=[1 3 5 8 10 12 15 17 19];
K=[2 4 6 7 9 11 13 14 16 18];
%find the strain variations
DU.DUxxc=acqEc.Uxx-repmat(e.pre.Uxx,length(acqEc.Uxx),1);
DU.DUyyc=acqEc.Uyy-repmat(e.pre.Uyy,length(acqEc.Uxx),1);
DU.DUxyc=acqEc.Uxy-repmat(e.anl.Uxyf,length(acqEc.Uxx),1);

%find the strain variations for each components
DU.DU1c(:,V)=(DU.DUxxc(:,V)+DU.DUyyc(:,V)+2*DU.DUxyc(:,V))/2;
DU.DU3c(:,V)=(DU.DUxxc(:,V)+DU.DUyyc(:,V)-2*DU.DUxyc(:,V))/2;
DU.DU2c(:,V)=DU.DUyyc(:,V);

DU.DU1c(:,K)=(DU.DUxxc(:,K)+DU.DUyyc(:,K)+2*DU.DUxyc(:,K))/2;
DU.DU3c(:,K)=(DU.DUxxc(:,K)+DU.DUyyc(:,K)-2*DU.DUxyc(:,K))/2;
DU.DU2c(:,K)=DU.DUyyc(:,K);

%find the strain variations
DU.DUxx=acqE.Uxx-repmat(e.pre.Uxx,length(acqE.Uxx),1);
DU.DUyy=acqE.Uyy-repmat(e.pre.Uyy,length(acqE.Uxx),1);
DU.DUxy=acqE.Uxy-repmat(e.anl.Uxyf,length(acqE.Uxx),1);

%find the strain variations for each components
DU.DU1(:,V)=(DU.DUxx(:,V)+DU.DUyy(:,V)+2*DU.DUxy(:,V))/2;
DU.DU3(:,V)=(DU.DUxx(:,V)+DU.DUyy(:,V)-2*DU.DUxy(:,V))/2;
DU.DU2(:,V)=DU.DUyy(:,V);

DU.DU1(:,K)=(DU.DUxx(:,K)+DU.DUyy(:,K)+2*DU.DUxy(:,K))/2;
DU.DU3(:,K)=(DU.DUxx(:,K)+DU.DUyy(:,K)-2*DU.DUxy(:,K))/2;
DU.DU2(:,K)=DU.DUyy(:,K);


DU.intVdtMat=acqEc.intVdtMat;
e.DU=DU;

%-----repair location by xcorr
% A=get_A_at_x(phE,acqE(1).x_sg,1);
%
%  for j=1:19 A.tv(:,j)=(A.t-front.tSg(j)).*front.vSg(j); end
% for j=1:19 acqE.tv(:,j)=(acqE.t-front.tSg(j))*front.vSg(j); end
% tv_lag=calc_xcorr_from_normlized_A(A.tv,A.lines,-12,7);

%--create t*Vfront axis for acqE and A
%A.tv=A.tv+repmat(tv_lag,length(A.tv),1);


figure;
plot(frontRaw.xv,frontRaw.v,'.-');
% hold all;
% plot(front.xv,front.v,'-');
% plot(frontSg.x,frontSg.v,'o-');
%plot(A.intVdt(2:end),diff(A.intVdt)./diff(A.t),'.-');
%plot(A.intVdt(A.intVdt_i1:A.intVdt_i2-1),diff(A.intVdt(A.intVdt_i1:A.intVdt_i2))./diff(A.t(A.intVdt_i1:A.intVdt_i2)),'.-');

function [xv,v]=calc_front_v(x,t)

if(x(1)<x(end))
    [~,index]=min(abs(x-50));
    [~,index2]=min(abs(x-90));
else
    [~,index]=min(abs(x-120));
    [~,index2]=min(abs(x-80));
end

vtmp=smooth(diff(x)./diff(t),3)';
v1=vtmp(1:index);
vtmp=smooth(diff(x)./diff(t),3)';
v2=vtmp(index+1:index2);
vtmp=smooth(diff(x)./diff(t),3)';
v3=vtmp(index2+1:end);
v=[v1 v2 v3];
xv=x(2:end)-0.5*(mean(diff(x)));
