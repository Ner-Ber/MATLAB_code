function anls_BiMaterial(e)

%e=anls_front_in_space(exper{6},7,30,190);
% sg_index1=6;
% sg_index2=14;
sg_index=[6:11 13 15];%
%sg_index=[6:11 13];%15 16
%sg_index=[(5:11)]; %5.5mm on 30mm
%sg_index=[3:5 9 13]; %7.5mm on 30mm
x_sg=e.acqE.x_sg(sg_index);

fig_index=20;
figure(fig_index);
fig_index=fig_index+1;
hold all;
plot(e.frontRaw.xv,e.frontRaw.v,'.-');
my_legend_add([e.acqE.Date '\' e.acqE.exp '\' num2str(e.acqE.event) ]);
plot(e.front.xv,e.front.v,'.-');
my_legend_add([e.acqE.Date '\' e.acqE.exp '\' num2str(e.acqE.event) ]);

figure(fig_index);
fig_index=fig_index+1;
hold all;
plot(x_sg,e.anl.UxxP(sg_index),'.-');
my_legend_add([e.acqE.Date '\' e.acqE.exp '\' num2str(e.acqE.event) ]);
plot(x_sg,e.anl.UyyP(sg_index),'.-');
my_legend_add([e.acqE.Date '\' e.acqE.exp '\' num2str(e.acqE.event) ]);

%---------------calc Af & Aovershoot automatically;
figure(fig_index);
fig_index=fig_index+1;
hold all;
plot(e.phECut.frontX,e.phECut.Af,'.-');
plot(e.phECut.frontX,e.phECut.AOverShoot,'.-');
Af=smooth(e.phECut.Af,5);
AOverShoot=smooth(e.phECut.AOverShoot,5);
plot(e.phECut.frontX,Af,'.-');
plot(e.phECut.frontX,AOverShoot,'.-');

for j=1:length(x_sg) [~,index(j)]=min(abs(e.phECut.frontX-x_sg(j))); end
Af=Af(index);
AOverShoot=AOverShoot(index);
plot(e.phECut.frontX(index),Af,'o-');
plot(e.phECut.frontX(index),AOverShoot,'o-');

% %-------calc Af & Aovershoot manually from snaphots;
figure(fig_index);
fig_index=fig_index+1;
hold off;
for j=1:length(x_sg) 
    
    [~,index]=min(abs(e.phECut.frontX-x_sg(j))); 
    
plot(e.phECut.xOffset(index,:), e.phECut.lines(index,:),'.-');  
my_legend_add(num2str(x_sg(j)));
plot(e.phECut.xOffsetSmt{index}, e.phECut.linesSmt{index},'.-');  
hold all;
xlim([-15 5]);
my_legend_add(num2str(x_sg(j)));
end
%[~,y]=my_ginput;
%Af=y(1,:);
%AOverShoot=y(2,:);


% %-----------calc Af & Aovershoot manually from same location;
% figure(fig_index);
% fig_index=fig_index+1;
% hold off;

 plot(e.A.intVdtMat(:,sg_index),subtruct_norm(e.A.lines(:,sg_index),1),'.-');
 xlim([-15 5]);
 my_legend_add(num2str(x_sg'));

%[~,y]=my_ginput;
%Af=y(1,:);
%AOverShoot=y(2,:);


%-------------Plot A vs Uxx
figure(fig_index);
fig_index=fig_index+1;
hold all;
plot(e.anl.UxxP(sg_index), 1-AOverShoot./Af,'.-');
my_legend_add([e.acqE.Date '\' e.acqE.exp '\' num2str(e.acqE.event) ]);
xlabel('Uxx');
ylabel('1-AOverShoot./Af');



