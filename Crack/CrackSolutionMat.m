function O=CrackSolutionMat(v)
%Uses CrackSolutionForh(v,h) to calculate stress\strain mesh 


%h_vec=[linspace(-7,-0.1,100)*1E-3, linspace(0.1,7,100)*1E-3];%[m]
%dx=abs(h_vec(1)-h_vec(2));
% for broberg

% sol=CrackSolution_SelfSimilarFreund(0.965,100e-3,-3.5e-6);
% y0=1e-6;
% ymax=sol.l/sol.v*sol.Cd;
h_vec=linspace(1e-5,30e-3,200);
% O.l=sol.l;
% O.Cs=sol.Cs;



for j=1:length(h_vec)
j/length(h_vec)
h=h_vec(j);
%sol=CrackSolution_SelfSimilarFreund(0.96,100e-3,-h);
sol=CrackSolutionForh(0.1,-1/2,h);
%sol=Crack_SupershearCohesiveForhNew(2200,h);
%O.vx(j,:)=sol.vx;
%O.vy(j,:)=sol.vy;
%sol=CrackSolutionForh(v,-1/2,h);
 
O.Sxy(j,:)=sol.Sxy;
O.Sxx(j,:)=sol.Sxx;
O.Syy(j,:)=sol.Syy;

O.Uxy(j,:)=sol.Uxy;
%planeStress
O.Uxx(j,:)=sol.Uxx;
O.Uyy(j,:)=sol.Uyy;

end


O.h=h_vec;
O.x=sol.x;
O.v=sol.v;

figure;
% subplot(1,2,1);mesh(O.x,O.h,-O.Uxx);view(0,90);colorbar;
% subplot(1,2,2);mesh(O.x,O.h,O.Uxy);view(0,90);colorbar;
x=O.x*1e3;
h=[-O.h(end:-1:1) O.h]*1e3;
Uxy=[O.Uxy(end:-1:1,:); O.Uxy]*1e3;

%mesh(x(1:2:end),h,Uxy(:,1:2:end));colorbar;
view([0,90])
imagesc(x(1:2:end),h,Uxy(:,1:2:end));
set(gca,'ydir','normal');
h = colorbar;
%title(h,'U_{xy}')

figure;
% subplot(1,2,1);mesh(O.x,O.h,-O.Uxx);view(0,90);colorbar;
% subplot(1,2,2);mesh(O.x,O.h,O.Uxy);view(0,90);colorbar;
x=O.x*1e3;
h=[-O.h(end:-1:1) O.h]*1e3;
Uxx=[O.Uxx(end:-1:1,:); -O.Uxx]*1e3;

%mesh(x(1:2:end),h,Uxx(:,1:2:end));colorbar;
view([0,90])
imagesc(x(1:2:end),h,Uxx(:,1:2:end));
set(gca,'ydir','normal');
h = colorbar;
