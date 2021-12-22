function O=CrackSolutionSmoothedSg(v,h)
%v should be fraction of cr. i.e v=0.9 

%-----vishay sg parameters EA-06-015RJ-120

sg.h0=h;%[m]
sg.alpha=pi/4;

sg.l=0.38E-3;%[m] 
sg.width=0.51E-3;%[m] %approxiamtly square

dx=1E-5;%[m] 
dh=dx;
h_vec=sg.h0+1.9*sg.l:-dh:sg.h0-1.9*sg.l;%[m] 

%---create temp field matrices
for j=1:length(h_vec)

h=h_vec(j);
sol=CrackSolutionForh(v,-1/2,h);

tmp.Uxy(j,:)=sol.Uxy;

 tmp.Uxx(j,:)=sol.Uxx;
 tmp.Uyy(j,:)=sol.Uyy;

% % planeStrain
% tmpPstrn.Uxx(j,:)=sol.UxxPstrn; 
% tmpPstrn.Uyy(j,:)=sol.UyyPstrn;

end
tmp.x=sol.x;
clear sol;

%-smooth
[O.Uxx O.Uxy O.Uyy]=smooth_sg(tmp,sg,h_vec);
% tmp.Uxx=tmpPstrn.Uxx;
% tmp.Uyy=tmpPstrn.Uyy;
% [O.UxxPstrn O.Uxy O.UyyPstrn]=smooth_sg(tmp,sg,h_vec);

O.x=tmp.x;
O.v=v;

function [Uxx Uxy Uyy]=smooth_sg(tmp,sg,h_vec)

%------upper gage
hTop=sg.h0+sg.l/2+sg.l;
hBottom=sg.h0+sg.l/2;

[~,index_hTop]=min(abs(h_vec-hTop));
[~,index_hBottom]=min(abs(h_vec-hBottom));

%---Smooth Uyy
dx=abs(tmp.x(1)-tmp.x(2));
dh=abs(h_vec(1)-h_vec(2));

smt_x=ceil(sg.width/dx);
Uyy=mean(tmp.Uyy(index_hTop:index_hBottom,:),1);%smooth for different h
Uyy=smooth(Uyy,smt_x)';

%---2 lower gages - calculate transverse parts of the rossete

hTop=sg.h0-sg.l/2;
hBottom=sg.h0-sg.l/2-sg.l*cos(sg.alpha);

[~,index_hTop]=min(abs(h_vec-hTop));
[~,index_hBottom]=min(abs(h_vec-hBottom));

sg_minus=(tmp.Uxx+tmp.Uyy+2*tmp.Uxy)/2; %two bottom legs of the rossete  
sg_minus=sg_minus(index_hTop:index_hBottom,:);
sg_plus=(tmp.Uxx+tmp.Uyy-2*tmp.Uxy)/2;
sg_plus=sg_plus(index_hTop:index_hBottom,:);

%---Smooth 

for j=2:length(sg_plus(:,1))
shift=floor(tan(sg.alpha)*dh*(j-1)/dx);
sg_plus(j,:)=circshift(sg_plus(j,:),[1,-shift]); %pointing to positive direction
sg_minus(j,:)=circshift(sg_minus(j,:),[1,shift]); %pointing to negative direction
end
sg_plus=mean(sg_plus,1);
sg_minus=mean(sg_minus,1);

smt_x=ceil(sg.width/cos(sg.alpha)/dx); %effective width in x direction is larger
sg_plus=smooth(sg_plus,smt_x)';
sg_minus=smooth(sg_minus,smt_x)';

%--calc Uxx,Uxy
Uxy=(sg_minus-sg_plus)/2;
Uxx=(sg_minus+sg_plus)-Uyy; %calculated transvers strain

