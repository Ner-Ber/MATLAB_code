function [a_vec,b_vec]=Crack_SupershearEqMotionCalcV_General(sol,x,tau_0)
% use [sol]=Crack_SupershearEqMotion3_General as an input
%[x] m, [tau]=tau_d

[Cd, Cs ,~, ~ ,~ ,~, mu, Gamma ,~,tau_d, Xc0]=CrackSolutionMaterialProperties;


tau_0=tau_0*tau_d;
a_tmp=linspace(55e-3,180e-3,100);

%---only above bmin
%bmin=sol.k*sqrt(2);
bmin=0.8;
g_vec=sol.g(sol.b>bmin);
G_my_vec=sol.G_my(sol.b>bmin);
b_vec=sol.b(sol.b>bmin);
% %---only below bmax
bmax=0.9;
g_vec=g_vec(b_vec<bmax);
G_my_vec=G_my_vec(b_vec<bmax);
b_vec=b_vec(b_vec<bmax);
%---subsample
g_vec=g_vec(1:2:end);
G_my_vec=G_my_vec(1:2:end);
b_vec=b_vec(1:2:end);

a_vec=g_vec*0;



parfor j=1:length(g_vec)
    
    g=g_vec(j);
    G_my=G_my_vec(j);
    Ks=@(a)calc_Ks(x,tau_0,g,a);
    
    fun_tmp=@(a) Gamma-(Ks(a)).^(1/g).*tau_d^2/mu/(tau_d).^(1./g).*G_my; %zero of this function determines a
    
    fun_tmp_vec=fun_tmp(a_tmp);
    [maxValue]=max(fun_tmp_vec);
    [minValue]=min(fun_tmp_vec);
    [~,index]=min(abs(fun_tmp(a_tmp)));
    
    if ((maxValue>0)&&(minValue<0))%zero exists
        
        a_vec(j)=fzero(fun_tmp,a_tmp(index));
    else
        a_vec(j)=NaN;
        
    end
    
end



function Ks=calc_Ks(x,tau_0,g,a)

%Ks=tau_0*a.^g/g;

tau_0=@(s) interp1(x,tau_0,s);
Ks=a*0;
for j=1:length(a)
    I=@(a,s) tau_0(a-s)./(s).^(1-g);
    Ks(j)= my_integralSingular(@(s)I(a(j),s),1-g,0,a(j));
end

