%  clear

n0=-0.4;
n=[n0-1,n0,n0+1];
 for j=1:3
    % sol{j}=CrackSolutionForh(0.965,j-1,3.5E-3);
     sol{j}=CrackSolutionSmoothedByMaskSg(0.96,n(j),3.5E-3);
 end

%H=100E-3;
H=50E-3;
%a=H^-0.5*1000*[-0.005,0*5,0*15];%n0=-0.25
a=H^-0.5*0.35*1000*[-0.025,3,1];%n0=-0.5
%a=H^-0.5*1000*[-0.05,10,50];%n0=-0.1
soltmp.x=sol{1}.x;
soltmp.Uxy=0*sol{1}.Uxy;%-a(1)*0.05E-3;
soltmp.Uxx=0*sol{1}.Uxx;
soltmp.Uyy=0*sol{1}.Uyy;

%---------separate plot
for j=1:3
    soltmp.Uxy=a(j)*H^(-sol{j}.n)*sol{j}.Uxy;
    soltmp.Uxx=a(j)*H^(-sol{j}.n)*sol{j}.Uxx;
    soltmp.Uyy=a(j)*H^(-sol{j}.n)*sol{j}.Uyy;
    
figN=2;
figure(figN);hold all;
plot(soltmp.x*1000,soltmp.Uxy,'-');
ch=get(gca,'Children');
c=get(ch,'color');
c=cell2mat(c);
figure(figN+1);hold all;
plot(soltmp.x*1000,soltmp.Uxx,'-','color',c(1,:));
figure(figN+2);hold all;
plot(soltmp.x*1000,soltmp.Uyy,'-','color',c(1,:));
    
end
%------------plot of the sum
for j=1:3
    soltmp.Uxy=soltmp.Uxy+a(j)*H^(-sol{j}.n)*sol{j}.Uxy;
    soltmp.Uxx=soltmp.Uxx+a(j)*H^(-sol{j}.n)*sol{j}.Uxx;
    soltmp.Uyy=soltmp.Uyy+a(j)*H^(-sol{j}.n)*sol{j}.Uyy;
end
% figN=2;
% figure(figN);hold all;
% plot(soltmp.x*1000,soltmp.Uxy,'.-');
% ch=get(gca,'Children');
% c=get(ch,'color');
% c=cell2mat(c);
% figure(figN+1);hold all;
% plot(soltmp.x*1000,soltmp.Uxx,'.-','color',c(1,:));
% figure(figN+2);hold all;
% plot(soltmp.x*1000,soltmp.Uyy,'.-','color',c(1,:));