function O=CrackSolutionSmoothedByMaskSg(v,h)
%Returns strain fields after smoothing by mask
%Units:  [v]=[Cr](should be afraction of Cr),
%[h]=[m], h should be positive
%n-the power for the LEFM solution. See definition in CrackSolution_n.m

sol=CrackSolutionForh(v,-1/2,h);
%sol=Crack_SupershearCohesiveForhNew(v,h);
%sol=Crack_SupershearCohesive(sym(1/7),h);
%sol=CrackSolution_SelfSimilarFreund(v,110e-3,-3.5e-3);

dx=sol.x(2)-sol.x(1);
sg.h0=h;%[m]
%dx=10E-6;%[m] Better choose as dx
dh=dx;
mask=CrackSolutionSgMask(dx);
l=length(mask.sgMinus(:,1));%should be odd
h_vec=sg.h0+(-(l-1)/2:1:(l-1)/2)*dh;
%h_vec=sg.h0+1.9*sg.l:-dh:sg.h0-1.9*sg.l;%[m]

%---create temp field matrices
for j=1:length(h_vec)
    %j/length(h_vec)
    h=h_vec(j);
    
    %-------Slip Weakening
    %     solTmp=CrackSolutionSW(v,h,2.5E-2);
    %     sol=solTmp.Sw;
    %     sol.x=solTmp.x;
    %     sol.Cr=solTmp.Cr;
    %
    %-------Slip pulse
    %     solTmp=CrackSolutionSP(v,2.5E-3,10,h);
    %     sol=solTmp.Sp;
    %     sol.x=solTmp.x;
    %     sol.Cr=solTmp.Cr;
    %-------LEFM
    sol=CrackSolutionForh(v,-1/2,h);
    %-----Broberg
    
    %sol=CrackSolution_SelfSimilarFreund(v,110e-3,-h);
%     %------SuperShear
%     j/length(h_vec)
%     sol=Crack_SupershearCohesiveForhNew(v,h);
% %    sol=Crack_Supershear(v,h);
%     
%      %sol=Crack_SupershearCohesive(sym(1/7),h);
     

tmp.Uxy(j,:)=sol.Uxy;
tmp.Uxx(j,:)=sol.Uxx;
tmp.Uyy(j,:)=sol.Uyy;
     
end
tmp.x=sol.x;

% O.v=sol.v;
% O.Cr=sol.Cr;
% O.Cd=sol.Cd;
% O.Cs=sol.Cs;
% O.n=sol.n;

O=rmfield(sol,{'Sxx','Syy','Sxy'});

%-smooth by using the mask
[O.sgMinus O.sgPlus O.Uxx O.Uxy O.Uyy O.x]=smooth_sg(tmp,mask);
O.y=sg.h0;
O.t=-O.x/O.v;

%--smooth to compinsate for finite time measurement
% smtT=1E-6*O.v;
% O.smtT=ceil(smtT/dx);
% 
% %-----------------calc Stress
% %---------Plane Stress (Szz=0)
% O.Sxx=E/(1-sigma^2)*(O.Uxx+sigma*O.Uyy); 
% O.Syy=E/(1-sigma^2)*(O.Uyy+sigma*O.Uxx);  
% O.Sxy=E/(1+sigma)*O.Uxy;
% note='Plane Stress';
% 
% %---------Plane Strain (Uzz=0)
% 
% % Sxx=-E/(1+sigma)/(1-2*sigma)*((1-sigma)*Uxx+sigma*Uyy); % "-" is added positive stress is compresion
% % Syy=-E/(1+sigma)/(1-2*sigma)*((1-sigma)*Uyy+sigma*Uxx); % "-" is added positive stress is compresion
% % Sxy=E/(1+sigma)*Uxy;
% % note='Plane Strain';
% % %----


function [sgMinus sgPlus Uxx Uxy Uyy x]=smooth_sg(tmp,mask)

%-calculate the strains in direction of the gages
sgMinusMat=(tmp.Uxx+tmp.Uyy+2*tmp.Uxy)/2; %two bottom legs of the rossete
sgPlusMat=(tmp.Uxx+tmp.Uyy-2*tmp.Uxy)/2;

%--pre allocation
x=zeros(length(tmp.Uyy(1,:))-mask.size,1);
Uyy=zeros(length(x),1);
sgMinus=zeros(length(x),1);
sgPlus=zeros(length(x),1);
x=zeros(length(x),1);

%-- smooth by using the mask

for j=1:length(x)
    Uyy(j)=sum(sum(tmp.Uyy(:,j:j-1+mask.size).*mask.sgUp))/mask.sgUpN;
    sgMinus(j)=sum(sum(sgMinusMat(:,j:j-1+mask.size).*mask.sgMinus))/mask.sgMinusN;
    sgPlus(j)=sum(sum(sgPlusMat(:,j:j-1+mask.size).*mask.sgPlus))/mask.sgPlusN;
    x(j)=tmp.x(j+(mask.size-1)/2);
end

%--calc Uxx,Uxy
Uxy=(sgMinus-sgPlus)/2;
Uxx=(sgMinus+sgPlus)-Uyy; %calculated transvers strain



