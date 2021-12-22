function [e]=anls_front_in_space(exper,event,Xstart,Xend,kind)

%Also works for rupture propagating in negative direction. For bi directional propagation use seperatly.
%kind- Rupture mode - Slow, Rayleigh, Super

exp_details=expDetailsRead(exper);

if (nargin<5)
    kind='Rayleigh';
end

if (strcmp(kind,'Slow'))
    smtPhE=11;
    smtAcq=7;%7
elseif(strcmp(kind,'Rayleigh'))
    smtPhE=5;
    smtAcq=1;
elseif(strcmp(kind,'Super'))
    smtPhE=3;
    smtAcq=1;
end

lines=3:6;

%-------Find tStart and tEnd
phE=phantomGetLines(exper,event,-6,'min',6,1000,smtPhE,'',lines);
A=get_A_at_x(phE,[Xstart, Xend],5);
index=find(subtruct_norm(A.lines(:,1),1)<0.95,1,'first');
tStart=A.t(index);
index=find(subtruct_norm(A.lines(:,2),1)<0.95,1,'first');
tEnd=A.t(index);

%------------Get the raw data

phE=phantomGetLines(exper,event,tStart-1,'min',tEnd+2,1000,smtPhE,'',lines);
phE.firstLine=mean(phE.lines(1:20,:),1);
acqE=acq132_event_get_data(exper,event,'start','end',smtAcq,'U1','U2','U3','x_sg','y_sg','sg_angle','host');

%acqE=acq132_event_get_dataVK(exper,event,'start','end',1,'Uxx','Uxy','Uyy','Sxx','Sxy','Syy','x_sg');
%acqE=mixVK_acq132_event_get_data(exper,event,'start','end',1,'Uxx','Uxy','Uyy','Sxx','Sxy','Syy','x_sg');

%---Fix for shear sensitivity

%%-Shift U1&U3

sub_sample_index=100;
dt=(acqE.t(2)-acqE.t(1))/sub_sample_index;%(msec)
if(strcmp(kind,'Super'))
    shift_t=1.6e-4;%0.3mm/2200m/s=0.14mus
else
    shift_t=2.5e-4;%0.3mm/1200m/s=0.25mus
end

d=ceil(shift_t/dt);

t_spline=(acqE.t(1):dt:acqE.t(end));
for j=1:length(acqE.x_sg)
    if (exp_details.sg_angle(j)==0) %Otherwise some thing more elaborated should be done.
        U1_spline=spline(acqE.t,acqE.U1(:,j),t_spline)';
        U2_spline=spline(acqE.t,acqE.U2(:,j),t_spline)';
        U3_spline=spline(acqE.t,acqE.U3(:,j),t_spline)';
        
        U1_spline=circshift(U1_spline,d);
        U3_spline=circshift(U3_spline,-d);
        
        acqE.U1(:,j)=U1_spline(1:sub_sample_index:end);
        acqE.U2(:,j)=U2_spline(1:sub_sample_index:end);
        acqE.U3(:,j)=U3_spline(1:sub_sample_index:end);
    end
end
disp(['shift=' num2str(shift_t)] )

acqE.gV=[0,0.1,0.95,-0.08];
%acqE.gV=[0.0,0.15,0.95,-0.08];
%acqE.gV=[0,0,1,0];
[acqE.U1,acqE.U2,acqE.U3]=calc_shear_sensitivity4(acqE.U1,acqE.U2,acqE.U3,acqE.gV);
disp('shear sensitivity was corrected' )
[~,~,~,acqE.Uxx,acqE.Uyy,acqE.Uxy,~]=calculate_stress_strain(acqE.U1,acqE.U2,acqE.U3,acqE.sg_angle);

%-----------------
A=get_A_at_x(phE,acqE(1).x_sg,1);

%----detect and smooth the front

[frontRaw front]=calc_frontXT_from_XAN(phE,tStart,tEnd,7,Xstart,Xend);

%-calc front v (set the smooth parameters at the function

[frontRaw.xv frontRaw.v]=calc_front_v(frontRaw.x,frontRaw.t,5); %smooth over 5mm

[front.xv front.v]=calc_front_v(front.x1,front.t,5);


%----------
%x=frontRaw.x; %another option is by front.x1
%t=frontRaw.t;

x=front.x1; %another option is by front.x1
t=front.t;

%--------convert Phe temporal to spatial measuremnts --- no need just use
%columns of phEcut.xOffset
% phE.intVdt=frontX_interp(x,t,0.9999999,phE.t);%0.999999
% [~,phE.intVdt_i1]=min(abs(phE.t-t(1)));
% [~,phE.intVdt_i2]=min(abs(phE.t-t(end) ));
%
% [~,index1]=min(abs(phE.x-x(1) ));
% [~,index2]=min(abs(phE.x-x(end) ));
% phE.intVdtMat=zeros(length(phE.t),length(phE.x));
%
% for j=index1:index2 phE.intVdtMat(:,j)=(-phE.intVdt+phE.x(j)); end



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

%------Find Af from time that was converted to space (columns of xOffset).

dx=5;
phECut.dx=dx;

phECut.xf=-25;
[~,index1]=min(abs(phECut.xOffset-(phECut.xf-dx)));
[~,index2]=min(abs(phECut.xOffset-(phECut.xf+dx)));
for j=1:length(index1)
    phECut.Af_from_t(j)=mean(phECut.lines(index2(j):index1(j),j),1);
end

phECut.xf2=-15;

[~,index1]=min(abs(phECut.xOffset-(phECut.xf2-dx)));
[~,index2]=min(abs(phECut.xOffset-(phECut.xf2+dx)));
for j=1:length(index1)
    phECut.Af2_from_t(j)=mean(phECut.lines(index2(j):index1(j),j),1);
end

phECut.xf3=-35;
[~,index1]=min(abs(phECut.xOffset-(phECut.xf3-dx)));
[~,index2]=min(abs(phECut.xOffset-(phECut.xf3+dx)));
for j=1:length(index1)
    phECut.Af3_from_t(j)=mean(phECut.lines(index2(j):index1(j),j),1);
end
clear index1 index2;

%---calc Xc

lines=phECut.lines-repmat(phECut.Af2_from_t,size(phECut.xOffset,1),1);
lines=lines./(1-repmat(phECut.Af2_from_t,size(phECut.xOffset,1),1));
[~,index1]= min(abs(lines-0.37));
for j=1:size(phECut.xOffset,2)
    phECut.Xc(j)=phECut.xOffset(index1(j),j);
end
clear index1;

%------Find Af from space .

[~,index1]=min(abs(phECut.xOffset'-(phECut.xf2-dx)));
[~,index2]=min(abs(phECut.xOffset'-(phECut.xf2+dx)));
for j=1:length(index1)
    phECut.Af2_from_x(j)=mean(phECut.lines(j,index1(j):index2(j)),2);
end

%%%---------Smooth A profiles after shifting with respect to xtip
% smt_index=2;%smt over 2*smt_index
% for j=1:length(phECut.xOffset(:,1) )
%
%     if (j<smt_index+1)
%         j1=0;
%     else
%         j1=smt_index;
%     end
%
%     if (j>length(phECut.xOffset(:,1) )- smt_index)
%         j2=0;
%     else
%         j2=smt_index;
%     end
%
%     [xTmp{j},ATmp{j}]=averageVec(phECut.xOffset(j-j1:j+j2,:)',phECut.lines(j-j1:j+j2,:)');
%
% end
%
% phECut.xOffsetSmt=xTmp;
% phECut.linesSmt=ATmp;


%%%---------Smooth A profiles after shifting with respect to xtip & find Af AOverShoot
% xf=front.xf;
% smt_index=2;%smt over 2*smt_index
% for j=1:length(phECut.xOffset(:,1) )
%
%     if (j<smt_index+1)
%         j1=0;
%     else
%         j1=smt_index;
%     end
%
%     if (j>length(phECut.xOffset(:,1) )- smt_index)
%         j2=0;
%     else
%         j2=smt_index;
%     end
%
%     [xTmp,ATmp]=averageVec(phECut.xOffset(j-j1:j+j2,:)',phECut.lines(j-j1:j+j2,:)');
%     [~,index1]=min(abs(xTmp)); %rutpure tip
%     [~,index2]=min(abs(xTmp--5));
%     [~,index3]=min(abs(xTmp--xf));
%     ATmp= smooth(ATmp,5);
%
%     phECut.xOffsetSmt{j}=xTmp;
%     phECut.linesSmt{j}=ATmp;
%
%     phECut.Af(j)=ATmp(index3);
%     phECut.AOverShoot(j)=min(ATmp(index2:index1));
%
%     % figure(30);
%     % hold all;
%     % plot(xTmp,ATmp,'.-');
%     % my_legend_add(num2str(phECut.frontX(j)));
%     % xlim([-15 5]);
% end

%---------

% for j=1:length(t)
%     [~, index_vec(j)]=min(abs(front.t-t(j)));
% end
% phECut.Af=front.Af( index_vec);

%----front at sg locations
SgNumber=length(acqE.x_sg);

for j=1:SgNumber
    [~,index2]=min(abs(frontRaw.xv-acqE.x_sg(j) - 5));
    [~,index1]=min(abs(frontRaw.xv-acqE.x_sg(j) +5));
    frontSg.v(j)=mean(frontRaw.v(index1:index2));
end

%---for bad quality "front" doesn't work
%for j=1:SgNumber [~,index(j)]=min(abs(front.xRaw-acqE.x_sg(j))); end
%frontSg.Af=front.Af(index);


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
A.tOffset=repmat(A.t,1,length(acqE.x_sg))-repmat(frontSg.t,length(A.t),1);

%%--------Single strain gages - tOffset Only on 7.5mm Block!
% clear index;
% for j=1:length(acqE.s.x_Uxx)
%     [~,index(j)]=min(abs(x-acqE.s.x_Uxx(j)));
% end
% front_t_Tmp=t(index);
% front_x_Tmp=acqE.s.x_Uxx;
% acqE.s.tOffset_Uxx=repmat(acqE.t,1,length(acqE.s.x_Uxx))-repmat(front_t_Tmp,length(acqE.t),1);
%
% clear index;
% for j=1:length(acqE.s.x_Uyy)
%     [~,index(j)]=min(abs(x-acqE.s.x_Uyy(j)));
% end
% front_t_Tmp=t(index);
% front_x_Tmp=acqE.s.x_Uyy;
% acqE.s.tOffset_Uyy=repmat(acqE.t,1,length(acqE.s.x_Uyy))-repmat(front_t_Tmp,length(acqE.t),1);


%--- Do better then constant rupture front velocity (Integrate Vdt)

acqE.intVdt=frontX_interp(x,t,0.9999999,acqE.t);%0.999999
[~,acqE.intVdt_i1]=min(abs(acqE.t-t(1)));
[~,acqE.intVdt_i2]=min(abs(acqE.t-t(end) ));

for j=1:SgNumber acqE.intVdtMat(:,j)=(-acqE.intVdt+frontSg.x(j)); end


A.intVdt=frontX_interp(x,t,0.9999999,A.t);
[~,A.intVdt_i1]=min(abs(A.t-tStart));
[~,A.intVdt_i2]=min(abs(A.t-tEnd));

%-----pre conditions
%e.pre=anls_pre_conditions(exper,event,-3,2:7,1); Doesn't use shear
%sensitivity correction
e.pre.Uxx=mean(acqE.Uxx(50:150,:),1);
e.pre.Uyy=mean(acqE.Uyy(50:150,:),1);
e.pre.Uxy=mean(acqE.Uxy(50:150,:),1);
% e.pre.x_sg=acqE.x_sg;
% acqTmp.Uxx=e.pre.Uxx;
% acqTmp.Uxy=e.pre.Uxy;
% acqTmp.Uyy=e.pre.Uyy;
% [e.pre.Syy,e.pre.Sxx,e.pre.Sxy]=calcStressFromStrain(acqTmp,'false');


for j=1:SgNumber A.intVdtMat(:,j)=-A.intVdt+frontSg.x(j); end
%A=structCutTime(A,tStart,tEnd);
%--- Cut the data for IntVdt
acqEc=structCutTime(acqE,acqE.t(acqE.intVdt_i1),acqE.t(acqE.intVdt_i2));
Ac=structCutTime(A,A.t(A.intVdt_i1),A.t(A.intVdt_i2));

%------------------------

if(~strcmp(kind,'Super'))

            % %--------- Analyze Rayleigh ruptures
    for j=1:length(acqE.x_sg)
        
    % --------------------find Uxyf & Uxxf
        dx=5;
        e.anl.xf=phECut.xf;
        [~,index1]=min(abs(acqE.intVdtMat(:,j)-(e.anl.xf+dx)));
        [~,index2]=min(abs(acqE.intVdtMat(:,j)-(e.anl.xf-dx)));
        %[~,index1]=min(abs(acqE.tv(:,j)--35));
        %[~,index2]=min(abs(acqE.tv(:,j)--45));
        e.anl.Uxyf(j)=mean(acqE.Uxy(index1:index2,j),1);
        e.anl.Uxxf(j)=-(mean(acqE.Uxx(index1:index2,j),1)-e.pre.Uxx(j));
        
    %-------------------     
        e.anl.xf2=phECut.xf2;
        [~,index1]=min(abs(acqE.intVdtMat(:,j)-(e.anl.xf2+dx)));
        [~,index2]=min(abs(acqE.intVdtMat(:,j)-(e.anl.xf2-dx)));
        e.anl.Uxxf2(j)=-(mean(acqE.Uxx(index1:index2,j),1)-e.pre.Uxx(j));
                     
     %-----------------------------        
        e.anl.xf3=phECut.xf3;
        [~,index1]=min(abs(acqE.intVdtMat(:,j)-(e.anl.xf3+dx)));
        [~,index2]=min(abs(acqE.intVdtMat(:,j)-(e.anl.xf3-dx)));
        e.anl.Uxyf3(j)=mean(acqE.Uxy(index1:index2,j),1);
        e.anl.Uxxf3(j)=-(mean(acqE.Uxx(index1:index2,j),1)-e.pre.Uxx(j));
        

      
        
        %------Near tip vicinity  
        [~,index1]=min(abs(acqE.intVdtMat(:,j)-25));
        [~,index2]=min(abs(acqE.intVdtMat(:,j)--25));
        
        
        % --find Uxy2    (rutpure tip)
        %     [e.anl.Uxy2(j) , indexTip]=min(acqE.Uxy(index1:index2,j));
        %     indexTip=index1+indexTip-1;
        %     e.anl.Uxy2t(j)=acqE.t(indexTip);
        %    % -Do something slightly better - find center of mass
        %     Uxy0Tmp=mean(acqE.Uxy(indexTip-200:indexTip-100,j));
        %     UxyTmp=abs(acqE.Uxy(:,j)-Uxy0Tmp).^2;
        %     e.anl.Uxy2t(j)=sum( UxyTmp(indexTip-1:indexTip+1).*acqE.t(indexTip-1:indexTip+1) ) / ( sum( UxyTmp(indexTip-1:indexTip+1 ) ) );
        
        %   ------find Uxy1 (Bump location)
        %     [e.anl.Uxy1(j) , indexBump]=max(acqE.Uxy(index1:indexTip,j));
        %     indexBump=index1+indexBump-1;
        %     e.anl.Uxy1t(j)=acqE.t(indexBump);
        %
        %-------find Uxy0   - first arrival of the bump
        %     Uxy0Tmp=mean(acqE.Uxy(indexBump-200:indexBump-100,j));
        %     UxyTmp=smooth(acqE.Uxy(indexBump-100:indexBump,j),5)-Uxy0Tmp;
        %     UxyTmp=UxyTmp(end:-1:1);
        %     index0=find(UxyTmp<0,1,'first');
        %     index0=indexBump-index0+1;
        %     e.anl.Uxy0(j)=acqE.Uxy(index0,j);
        %     e.anl.Uxy0t(j)=acqE.t(index0);
        
        
        %%%------find pre close to rupture arrival
        % [~,index1]=min(abs(acqE.intVdtMat(:,j)-60));
        % e.pre.Uxx(j)=mean(acqE.Uxx(index1-20:index1,j));
        % e.pre.Uyy(j)=mean(acqE.Uyy(index1-20:index1,j));
        % e.pre.Uxy(j)=mean(acqE.Uxy(index1-20:index1,j));
        
        %----- find Peak
        [e.anl.UxxP(j) , indexP]= min(acqE.Uxx(index1:index2,j)-e.pre.Uxx(j));
        e.anl.UxxP(j)=-e.anl.UxxP(j);
        indexP=index1+indexP-1;
        e.anl.UxxPt(j)=acqE.t(indexP);
        e.anl.UxxPx(j)=acqE.intVdtMat(indexP,j);
        
        [e.anl.UyyP(j) , indexP]=max((acqE.Uyy(index1:index2,j)-e.pre.Uyy(j) ));
        indexP=index1+indexP-1;
        e.anl.UyyPt(j)=acqE.t(indexP);
        e.anl.UyyPx(j)=acqE.intVdtMat(indexP,j);
        
    end
    
else
    
    % % % %%%%---Analyze SuperShear
    for j=1:length(acqE.x_sg)
        
        [~,index1]=min(abs(acqE.tv(:,j)-25));
        [~,index2]=min(abs(acqE.tv(:,j)--25));
        
        %------find pre close to rupture arrival
        %         e.pre.Uxx(j)=mean(acqE.Uxx(index1-20:index1,j));
        %         e.pre.Uyy(j)=mean(acqE.Uyy(index1-20:index1,j));
        %         e.pre.Uxy(j)=mean(acqE.Uxy(index1-20:index1,j));
        
        %----- find Peak
        [e.anl.UxxP(j) , indexP]=max(abs(acqE.Uxx(index1:index2,j)-e.pre.Uxx(j)));
        indexP=index1+indexP-1;
        e.anl.UxxPt(j)=acqE.t(indexP);
        
        [e.anl.UyyP(j) , indexP]=max(abs(acqE.Uyy(index1:index2,j)-e.pre.Uyy(j))  );
        indexP=index1+indexP-1;
        e.anl.UyyPt(j)=acqE.t(indexP);
        
        %--------for mach cone calculations
        [e.anl.Uyy1(j) , indexP]=min( acqE.Uyy(index1:index2,j)-e.pre.Uyy(j)  );
        indexP=index1+indexP-1;
        e.anl.Uyy1t(j)=acqE.t(indexP);
        
        %--------find Uxyf
        [~,index1]=min(abs(acqE.tv(:,j)--30));
        [~,index2]=min(abs(acqE.tv(:,j)--40));
        %[~,index1]=min(abs(acqEc.intVdtMat(:,j)--45));
        %[~,index2]=min(abs(acqEc.intVdtMat(:,j)--55));
        %[~,index1]=min(abs(acqEc.intVdtMat(:,j)--20));
        % [~,index2]=min(abs(acqEc.intVdtMat(:,j)--40));
        e.anl.Uxyf(j)=mean(acqE.Uxy(index1:index2,j),1);
        e.anl.Uxxf(j)=-(mean(acqE.Uxx(index1:index2,j),1)-e.pre.Uxx(j));

    end
end

e.frontRaw=frontRaw;
e.frontSg=frontSg;
e.front=front;
e.acqE=acqE;
e.acqE.smtAcq=smtAcq;
e.acqEc=acqEc;
e.phE=phE;
e.phECut=phECut;
e.A=A;
e.Ac=Ac;



%-----repair location by xcorr
% A=get_A_at_x(phE,acqE(1).x_sg,1);
%
%  for j=1:19 A.tv(:,j)=(A.t-front.tSg(j)).*front.vSg(j); end
% for j=1:19 acqE.tv(:,j)=(acqE.t-front.tSg(j))*front.vSg(j); end
% tv_lag=calc_xcorr_from_normlized_A(A.tv,A.lines,-12,7);

%--create t*Vfront axis for acqE and A
%A.tv=A.tv+repmat(tv_lag,length(A.tv),1);

% fig_num=20;
% figure(fig_num);
% ax=findobj(fig_num,'Type','axes');%if fig_num is axes handle the command returns ax=fig_num
% axe=findobj(ax,'Type','line');
% DisplayName=get(axe,'DisplayName');
% DisplayName= [DisplayName(end:-1:1);num2str(e.phE.EventNum)];

%fig=my_get_axis;
%fig.DisplayName{end+1}=num2str(e.phE.EventNum);


% plot(frontRaw.xv,frontRaw.v,'.-');
% hold all;
%legend(DisplayName);legend off;

% hold all;
% plot(front.xv,front.v,'-');
% plot(frontSg.x,frontSg.v,'o-');
%plot(A.intVdt(2:end),diff(A.intVdt)./diff(A.t),'.-');
%plot(A.intVdt(A.intVdt_i1:A.intVdt_i2-1),diff(A.intVdt(A.intVdt_i1:A.intVdt_i2))./diff(A.t(A.intVdt_i1:A.intVdt_i2)),'.-');


function [xv,v]=calc_front_v(x,t,smtX)

xv=x(2:end)-0.5*(mean(diff(x)));

% if(x(1)<x(end))
%     [~,index]=min(abs(x-50));
%     [~,index2]=min(abs(x-90));
% else
%     [~,index]=min(abs(x-120));
%     [~,index2]=min(abs(x-80));
% end
%v=smooth(x(2:end),diff(x)./diff(t),5)';

v=diff(x)./diff(t);
v=my_smooth2(xv,v,smtX);



%vtmp=smooth(diff(x)./diff(t),3)';
% v1=vtmp(1:index);
% vtmp=smooth(diff(x)./diff(t),3)';
% v2=vtmp(index+1:index2);
% vtmp=smooth(diff(x)./diff(t),3)';
% v3=vtmp(index2+1:end);
% v=[v1 v2 v3];