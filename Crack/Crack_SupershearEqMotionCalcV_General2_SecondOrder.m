function [a_vec,b_vec1,b_vec2]=Crack_SupershearEqMotionCalcV_General2_SecondOrder(sol,x,tau_0,aZO,bZO)
%An attempt to take velocity variation in calculations of Ks based on zero
%order Cf profile.
%a-crack length, b-cf/cd
%[x] m, [tau]=tau_d
%aZO, gZO - zero order a and g
%solves for equation of motion for homogeneous material parameters for a given spatial profile of tau_0.

%[Cd, Cs ,~, ~ ,~ ,~, mu, Gamma ,~,tau_d, Xc0]=CrackSolutionMaterialProperties;
%a0=Gamma/tau_d^2*mu;

tau_d=1;%from here on every thing in normalized units
a0=1;
tau_0=tau_0*tau_d;
a_vec=linspace(1,max(x),150)*a0;%in uni0.0ts of a0

%-----calc gZO
aZO=[0 aZO];
bZO=[bZO(1) bZO];
k=sol.k;
ap=@(b)(1-b.^2).^0.5;
bs=@(b)((b/k).^2-1).^0.5;
g=@(b) 1/pi*atan(4*ap(b).*bs(b)./(bs(b).^2-1).^2);%g stands for gamma. Note that eq. 34 (book eq6.9.111) is for 1/b
gZO=g(bZO);
clear ap bs g;

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


for j=25:length(a_vec)-100
    a=a_vec(j);
    
   Ks=calc_Ks(x,tau_0,a,aZO,gZO);
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


a_vec=a_vec/a0;

function Ks=calc_Ks(x,tau_0,a,aZO,gZO)

%For homogeneous case Ks should be Ks=tau_0*a.^g. Therefore N_factor is
%used.

tau_0=@(s) interp1(x,tau_0,s);

%gZO=@(s) interp1(aZO,gZO,s*0+a);%Check - to get original result
gZO=@(s) interp1(aZO,gZO,s);

    %---- Semi-infinite crack
    I=@(s) tau_0(a-s)./(s).^(1-gZO(a-s));
    N_factor=1./gZO(a);
    
    %---- Bilateral crack
   % I=@(s) tau_0(a-s).*(a./(s.*(2*a-s))).^(1-gZO(a-s));
   % N_factor=pi^0.5*gamma(gZO(a))/(2*gamma(0.5+gZO(a)));
    
     Ks=my_integralSingular(I,1-gZO(a),0,a)/N_factor;


function b_zero=ZeroFromLinearInterp(b,fun_tmp,index)
%linear interp y=alpha*x+beta
            alpha=(fun_tmp(index+1)-fun_tmp(index-1))/(b(index+1)-b(index-1));
            beta=fun_tmp(index)-alpha*b(index);
            b_zero=-beta/alpha;