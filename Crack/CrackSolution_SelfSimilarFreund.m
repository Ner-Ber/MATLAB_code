function sol=CrackSolution_SelfSimilarFreund(v,l,y)

%Symmetric expansion of a shear crack -Freund 6.3.3 p330
% %v=[Cr], l=vt [m] , y=[m]
%use for x>0,y>0. if y<0 Sxy change sign-I don't understand why- probably some
%branch cut that I haven't taken into account
%-------------

[Cd Cs Cr ,~, ~, ~, mu G]=CrackSolutionMaterialProperties;

k=Cs/Cd;%Broberg p.330
v=v*Cr;
% -----------------calc beta
R=@(z) (Cs^-2 - 2*z.^2).^2 + 4*z.^2 .* (Cd^-2 -z.^2).^0.5.* (Cs^-2-z.^2).^0.5; %just after eq. 6.3.42
D=@(v) -v^4*real(R(1/v)); %Due to numerical error R has an imaginary part
%------First calc K
alpha_s=(1-(v/Cs).^2).^0.5;
%----- Following Broberg p.334,336
A=2*(1-k^2)*alpha_s.*v.^2./D(v)/Cs^2;
K=(G*mu*4*(1-k^2)./A).^0.5;


I_II_integrand=@(z) R(1i*z) .* (v^-2 + z.^2).^-1.5 .* (Cs^-2 + z.^2).^-0.5;
I_II= (Cs^2/v * integral( I_II_integrand, 0, Inf ) )^-1;

%---------take singular K and l then calculate tau_0
%for given crack length (l) and crack velocity, tau_0 is calculated to correspond to the intensity
%factor calculated before

tau_0= -K / ( I_II *  R(1/v) / (Cs^-2* v^-1 *(v^-2 - Cs^-2).^0.5) *(pi*l)^0.5 );  %eq6.3.67(Mistake?)
tau_0=real(tau_0);%Due to numerical error beta has an imaginary part

%---------take singular K and tau_0 then calculate l
% tau_0=0.2e-3*2*mu;%[Pa] %typical strain drop 0.1e-3
% l=( -K / ( I_II *  R(1/v) / (Cs^-2* v^-1 *(v^-2 - Cs^-2).^0.5) *tau_0*(pi)^0.5 ))^2;
% l=real(l);%Due to numerical error beta has an imaginary part
%---------

beta=real(tau_0/mu *I_II);%Due to numerical error beta has an imaginary part
%The integrals are written in normlized units : z->z*v, therfore a,b \
%are also normlized normlized units. x_i=x_i/(t*v)*)

z= @(x,y,v) x./(x.^2 + y.^2) +heaviside(1 - (x.^2 + y.^2) *(v)^2)* 1i*y./(x.^2 + y.^2).* (1 - (x.^2 + y.^2) *(v)^2).^0.5 -heaviside(-1 +(x.^2 + y.^2) *(v)^2)*1e-20i;% "eq.6.3.23";heaviside affects only for Gxt.
%z= @(x,y,v) x./(x.^2 + y.^2) + 1i*y./(x.^2 + y.^2).* (1 - (x.^2 + y.^2) *(v)^2).^0.5;% "eq.6.3.23";

%calculate the potentials. For F  v/Cd should be used,  for G use v/Cs
Fxy_1=@(z,v)  z .*(1 - z.^2).^0.5 .*(-v^2 + z.^2) ...
    - 2* v^2 *(-1 + z.^2).*(1 - z.^2/v^2).^0.5 .* ellipticE(asin(z), 1/v^2) ...
    +  v^2* (-1 + z.^2).* (1 - z.^2/v^2).^0.5.* ellipticF( asin(z), 1/v^2) ;
Fxy=@(z,v) 2i*Fxy_1(z,v) ./ ( (v^2 - z.^2).^0.5.*(-1 + z.^2));

Fxx=@(z,v) -2*(2 - z^2)/(1 - z.^2)^0.5;
Fyy=@(z,v) -2*(-2 + v^2 + z.^2)./(1 - z^2)^0.5;
Gxy= @(z,v) (4 - v^2 - 2*z.^2)./(1 - z.^2).^0.5;

Gxx_1=@(z,v) (-2 + v^2)*z.* (1 - z.^2).^0.5.*(-v^2 + z.^2) ...
    - v^2* (-4 + 3 *v^2)* (-1 + z.^2).*(1 - z.^2/v^2).^0.5.*ellipticE(asin(z), 1/v^2)...
    +2* v^2* (-1 + v^2)* (-1 + z.^2) .*(1 - z.^2/v^2).^0.5.*ellipticF(asin(z), 1/v^2);
Gxx=@(z,v) 1i*Gxx_1(z,v) ./ ( (v^2-1)*(v^2 - z.^2).^0.5.* (-1 + z.^2));

Gyy_1=@(z,v) -v^2 *(-4 + v^2) *(-1 + z.^2).*(1 - z.^2/v^2).^0.5.*ellipticE(asin(z), 1/v^2) ...
    + (-2 + v^2)* ( z.*(1 - z.^2).^0.5 .* ( -v^2 + z.^2) + v^2* (-1 + z.^2).* (1 - z.^2/v^2).^0.5.*ellipticF(asin(z), 1/v^2));
Gyy=@(z,v) 1i*Gyy_1(z,v)./((v^2 - z.^2).^0.5.*(-1 + z.^2));

Fxt=@(z,v)  2*(z./(1 - z.^2).^0.5 - asin(z));
Fyt=@(z,v) -2i* ((v^2 - z.^2).^0.5./(1 - z.^2).^0.5 - log( (1 - z.^2).^0.5 +(v^2 - z.^2).^0.5 ) );
Gyt=@(z,v) ((-2 + v^2)*z)./(1 - z.^2).^0.5 + 2*asin(z);

Gxt_1=@(z,v) - 1i*(  (-2 + v^2) *(v^2 - z.^2) + 2* (-1 + v^2)*(v^2 - z.^2).^0.5.* (-1 + z.^2).^0.5 .* atan(  (-1 + z.^2).^0.5./(v^2 - z.^2).^0.5) ) ;
Gxt=@(z,v) Gxt_1(z,v)./ (  (-1 + v^2) *(1 - z.^2).^0.5.* (v^2 - z.^2).^0.5 );

%calculate stress
%here real values of v should be used. x , y stand for  x/l y/l  (l=vt).

SxyTmp =@(x,y,v) beta*mu* (v/Cs)^-2 *real(2* Fxy(z(x, y, v/Cd), v/Cd) + Gyy(z(x, y, v/Cs), v/Cs) -Gxx(z(x, y, v/Cs), v/Cs));%6.3.21
Sxy =@(x,y,v)  tau_0 +SxyTmp(x, y, v);
Sxx=@(x, y, v)  beta*mu*(v/Cs)^-2* imag(  (Cd/Cs)^2* Fxx(z(x, y, v/Cd),v/Cd) + ((Cd/Cs)^2 - 2)* Fyy(z(x, y, v/Cd),v/Cd) + 2 *Gxy(z(x, y, v/Cs), v/Cs)  );%(*6.3.21*)
Syy =@ (x, y,v)  beta*mu* (v/Cs)^-2 *imag( ((Cd/Cs)^2 - 2)* Fxx(z(x, y, v/Cd), v/Cd) + (Cd/Cs)^2 *Fyy(z(x, y, v/Cd), v/Cd) - 2* Gxy(z(x, y, v/Cs), v/Cs) );%(*6.3.21*)

vx=@(x, y, v)  v*(v/Cs)^-2*beta*imag( Fxt(z(x, y, v/Cd), v/Cd) + Gyt(z(x, y, v/Cs), v/Cs) );

vyTmp=@(x, y, v) v*(v/Cs)^-2*beta*real(Fyt(z(x, y, v/Cd), v/Cd) - Gxt(z(x, y, v/Cs), v/Cs) );
vy=@(x, y, v)  -vyTmp(Cd/v,0,v)+vyTmp(x,y,v);


% %---------plot versus space
% if (l>50e-3)
%  %x=l-10e-3:10^-4:l/v*Cs+10e-3;
%   x=0.8*l/v*Cs:10^-4:l/v*Cd;
% %x=l-60e-3:1e-4:l+60e-3;
% else
%    x=l-0.5*l:10^-5:l/v*Cs+10e-3;
% 
% end
% 
% 
% warning off;
% for j=1:length(x)
% sol.Sxy(j)=Sxy(x(j)/l,y/l,v);
% sol.Sxx(j)=Sxx(x(j)/l,y/l,v);
% sol.Syy(j)=Syy(x(j)/l,y/l,v);
% 
% % sol.Sxy(j)=0;
% % sol.Sxx(j)=0;
% % sol.Syy(j)=0;
% 
%  sol.vx(j)=vx(x(j)/l,y/l,v);
%  sol.vy=vy(x/l,y/l,v);
%  end
% % 
%  warning on;
%  sol.x=x-l;
%  sol.smtT=ceil(1E-6/(abs(mean(diff(sol.x)))/v));


% % %---------plot polar
% r=1/v*Cs;
% % %r=1e-4:1e-3:1/v*Cd;
% %theta=(0.02:0.1:89.9)/180*pi;
%  theta=0.0001;
% 
% 
% warning off;
% 
% for j=1:length(theta)
%     for m=1:length(r)
%         sol.Sxy(m,j)=Sxy(r(m)*cos(theta(j)),r(m)*sin(theta(j)),v);
%         sol.Sxx(m,j)=Sxx(r(m)*cos(theta(j)),r(m)*sin(theta(j)),v);
%         sol.Syy(m,j)=Syy(r(m)*cos(theta(j)),r(m)*sin(theta(j)),v);
%         sol.vx(m,j)=vx(r(m)*cos(theta(j)),r(m)*sin(theta(j)),v);
%         sol.vy(m,j)=vy(r(m)*cos(theta(j)),r(m)*sin(theta(j)),v);
%     end
% end
% warning on;
% sol.r=r;
% sol.theta=theta;


%--------plot versus time
if (l>50e-3)
x=l-50e-3:0.5*10^-4:l+50e-3;
t=x/v;
else
    x=10e-3:0.5*10^-4:l+50e-3;
t=x/v;

end
t=t(end:-1:1);

warning off;
for j=1:length(x)
sol.Sxy(j)=Sxy(l/v/t(j),y/v/t(j),v);
sol.Sxx(j)=Sxx(l/v/t(j),y/v/t(j),v);
sol.Syy(j)=Syy(l/v/t(j),y/v/t(j),v);
sol.vx(j)=vx(x(j)/l,y/l,v);
sol.vy(j)=vy(x(j)/l,y/l,v);
end
warning on;
sol.t=t-l/v;
sol.smtT=ceil(1E-6 /abs (sol.t(2)-sol.t(1)) );
sol.x=-sol.t*v;

[sol.Uxx,sol.Uyy,sol.Uxy]=calcStrainFromStress(sol,mu,k);
sol.l=l;
sol.y=y;
sol.v=v;
sol.Cr=Cr;
sol.Cs=Cs;
sol.Cd=Cd;
sol.G=G;
sol.Uxy0=tau_0/2/mu;
sol.Sxy0=tau_0;
sol.mu=mu;

function [Uxx Uyy Uxy]=calcStrainFromStress(sol,mu,k)

%-------Works for both boundary conditions with appropriate k

Uxy=1/(2*mu)*sol.Sxy;
Uxx=1/(2*mu)/(2*(1-k^2))*(sol.Sxx-(1-2*k^2)*sol.Syy);
Uyy=1/(2*mu)/(2*(1-k^2))*(sol.Syy-(1-2*k^2)*sol.Sxx);


