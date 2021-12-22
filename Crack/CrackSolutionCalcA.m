function A=CrackSolutionCalcA(v,mode)
%[v]=m/s

[Cd Cs Cr nu ro E mu G PlaneStrain]=CrackSolutionMaterialProperties;

k=Cs/Cd;%Broberg p.330
alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;

if (nargin==1)
    mode=2;
end

if(mode==1)
    A=2*(1-k^2)*alpha_d.*v.^2./D/Cs^2; %ModeI
else
    A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2; %ModeII
end



