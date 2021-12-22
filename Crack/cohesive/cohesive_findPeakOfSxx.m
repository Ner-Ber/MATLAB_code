function cohesive_findPeakOfSxx()


[Cd, Cs, Cr, nu , ~, ~, mu, GammaDefault, PlaneStrain, tau_p_def]=CrackSolutionMaterialProperties(PlaneFlag);

k=Cs/Cd;%Broberg p.330
v=v*Cr;

%------- General LEFM functions
alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;

A0=-1i/mu/(2*pi)^0.5*2*alpha_s/D*K;%Eq.15


%------calc F
z=z_d;
syms x z
f = (-x).^0.5.*tau(x)./(x-(z+3.5e-3i));
% I=tau_p*calcI(f,z/L,ILowerBoundery)*L^0.5;
solz = f(-(z+3.5e-3i).^(-1.5).*(A0/2+(4*alpha_s./(1*pi*1i*mu*D)).*(1+(z+3.5e-3i)).*int(f,x,-L,0))==0,z)
F = A0*z.^-0.5+(2*pi*1i)^-1*4*alpha_s/mu/D*z.^-0.5.*I;%Eq.13





end

%% external functions
function I=calcI(f,z,ILowerBoundery)

for j=1:length(z)
    f_tmp=@(x) f(x,z(j));
    I(j)=integral(f_tmp,ILowerBoundery,0);
end


end