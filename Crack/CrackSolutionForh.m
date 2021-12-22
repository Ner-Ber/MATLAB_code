function sol=CrackSolutionForh(v,n,h,x)
%For constant y=h,range of x and v the function finds r,theta and calls CrackSolution
%Units:  [v]=[Cr](should be afraction of Cr),[h]=[m]
% n-the power for the LEFm solution. See definition in CrackSolution_n.m


%-------calc C_Rayleigh
% tmp=CrackSolution_n(900:1370,0,0,1);
% [~,index]=min(abs(tmp.D));
% Cr=tmp.v(index);
%-------------

if (nargin<4)
    dx=1E-5;
    x=-50E-3:dx:50E-3;%[m] distance from the crack tip ,theta=0;
    %dx=1E-4;
    %x=-30E-3:dx:30E-3;%[m] distance from the crack tip ,theta=0;
end

%------find theta
theta=atan(h./x);
if (h>0)
    %theta should be 0<theta<pi;
    index=find(theta<0);
    theta(index)=theta(index)+pi;
else
    %theta should be -pi<theta<0;
    index=find(theta>0);
    theta(index)=theta(index)-pi;
end

%---------find r
r=(h.^2+x.^2).^0.5;

%-----calc LEFM solution
sol=CrackSolution_n(v,n,theta,r);
sol.h=h;
sol.x=x;%[mm]
sol.t=-sol.x/sol.v;
sol.smtT=ceil(1E-6/(abs(mean(diff(sol.x)))/sol.v));
sol.Xc=0;
sol.tau_p=NaN;
sol.y=h;






