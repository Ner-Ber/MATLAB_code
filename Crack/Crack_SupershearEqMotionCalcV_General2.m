function [a_vec1,b_vec1,a_vec2,b_vec2]=Crack_SupershearEqMotionCalcV_General2(sol,x,tau_0)
%a-crack length, b-cf/cd
%[x] a0, [tau]=tau_d
%solves for equation of motion for homogeneous material parameters for a given spatial profile of tau_0.

%[Cd, Cs ,~, ~ ,~ ,~, mu, Gamma ,~,tau_d, Xc0]=CrackSolutionMaterialProperties;
%a0=Gamma/tau_d^2*mu;

tau_d=1;%from here on every thing in normalized units
a0=1;
tau_0=tau_0*tau_d;
a_vec=linspace(1,max(x),150)*a0;%in uni0.0ts of a0

%---only above bmin
%bmin=sol.k*sqrt(2);
bmin=0.5;
g=sol.g(sol.b>bmin);
G_my=sol.G_my(sol.b>bmin);
b=sol.b(sol.b>bmin);
%%---only below bmax
bmax=0.95;
g=g(b<bmax);
G_my=G_my(b<bmax);
b=b(b<bmax);

b_vec1=a_vec*0+NaN;
b_vec2=a_vec*0+NaN;


parfor j=1:length(a_vec)
   
    a=a_vec(j);
    %a=390;
    %g=1/2;
    Ks=calc_Ks(x,tau_0,a,g);%returns Ks for single value of a and vector of g
    fun_tmp=a0-(Ks/tau_d).^(1./g).*G_my; %zero of this function determines a
    
    [maxValue]=max(fun_tmp);
    [minValue]=min(fun_tmp);
    
    if ((maxValue>0)&&(minValue<0))%At least on zero exists
        
        [~,index]=findpeaks(-abs(fun_tmp));
        
        if (size(index,2)==1)
            
            b_vec1(j)=ZeroFromLinearInterp(b,fun_tmp,index);
                     
        elseif(size(index,2)==2)%If two branches exist
         
            b_vec1(j)=ZeroFromLinearInterp(b,fun_tmp,index(end));
            b_vec2(j)=ZeroFromLinearInterp(b,fun_tmp,index(1));
            
        end
        
    end
end


a_vec1=a_vec(~isnan(b_vec1))/a0;
b_vec1=b_vec1(~isnan(b_vec1));

a_vec2=a_vec(~isnan(b_vec2))/a0;
b_vec2=b_vec2(~isnan(b_vec2));

function Ks=calc_Ks(x,tau_0,a,g)

% For homogeneous case Ks should be Ks=tau_0*a.^g. Therefore N_factor is
% used.

tau_0=@(s) interp1(x,tau_0,s);
Ks=a*0;
for j=1:length(g)
    %---- Semi-infinite crack
   
   % I=@(s) tau_0(a-s)./(s).^(1-g(j));
   % N_factor=1./g(j);
    
    %---- Bilateral crack
     I=@(s) tau_0(a-s).*(a./(s.*(2*a-s))).^(1-g(j));
     N_factor=pi^0.5*gamma(g(j))/(2*gamma(0.5+g(j)));
    
     Ks(j)= my_integralSingular(I,1-g(j),0,a)/N_factor;
end

function b_zero=ZeroFromLinearInterp(b,fun_tmp,index)
%linear interp y=alpha*x+beta
            alpha=(fun_tmp(index+1)-fun_tmp(index-1))/(b(index+1)-b(index-1));
            beta=fun_tmp(index)-alpha*b(index);
            b_zero=-beta/alpha;