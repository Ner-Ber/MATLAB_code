function anls_EquationOfMotion_Plot2016_05_27_12_12_30(anl)


[Cd Cs Cr, ~ , ~ ,E mu Gamma, ~ ,~]=CrackSolutionMaterialProperties;

%Sligthly correct Cf
for j=1:length(anl.frontRaw)
    anl.frontRaw{j}.v=anl.frontRaw{j}.v*1.02;
end

index_sg=(1:length(anl.y3_5.x_sg)); %choose strain gages that will affect Uxy
index_sg_Plot=(1:length(anl.y3_5.x_sg)); %choose strain gages for plot
%index_sg_Plot=(4:length(anl.y3_5.x_sg)); %choose strain gages for plot
%index_sg=(5:length(anl.y3_5.x_sg));

%-----Define region to analyze

xStart=50;
xEnd=anl.y3_5.x_sg(index_sg(end));

%---------Define spatial variation of Gamma
anl.y3_5.Gamma=[anl.y3_5.x_sg]*0;
%---For 200mm
anl.y3_5.Gamma(anl.y3_5.x_sg>100)=Gamma;
anl.y3_5.Gamma(anl.y3_5.x_sg<100)=Gamma/1.67;
%G(anl.y3_5.x_sg(index_sg)<100)=Gamma;
anl.y3_5.Gamma(anl.y3_5.x_sg>137)=Gamma;
%---For 150mm
%G=G+Gamma;
%     G(anl.y3_5.x_sg>80)=Gamma;
%     G(anl.y3_5.x_sg<80)=Gamma/1.5;
%      G(anl.y3_5.x_sg<50)=Gamma/2;
%

%-------Correct Uxy0_Offset
sol=CrackSolutionForh(0.01,-1/2,3.5e-3);
[~, index]=min(abs(sol.x*1e3+40));
Uxy0_Offset=sol.Uxy(index)/Gamma^0.5*1e3;%First renormalize by Gamma
Uxy0_Offset=Uxy0_Offset.*anl.y3_5.Gamma.^0.5;%different locations might have different Gamma
anl.y3_5.Uxy0_f=anl.y3_5.Uxy0_f+repmat(Uxy0_Offset,length(anl.e),1);

% %-------calculate stuff
anl.y3_5.Cf=zeros(length(anl.e),length(anl.y3_5.x_sg));
anl.y3_5.Gs=zeros(length(anl.e),length(anl.y3_5.x_sg));

for j=1:length(anl.e)
    %for j=7:7
    %x=[0;19; anl.y3_5.x_sg(index_sg)']*1e-3;
    x=[0;35; anl.y3_5.x_sg(index_sg)']*1e-3;
    %-------Define pre stress
    %Uxy0=[0; anl.y3_5.Uxy0_f(j,index_sg)']*1e-3;
    Uxy0=[0;mean(anl.y3_5.Uxy0_f(j,index_sg(1:1))); anl.y3_5.Uxy0_f(j,index_sg)']*1e-3;
    %Uxy0=[mean(anl.y3_5.Uxy0_f(j,index_sg(1:3)));mean(anl.y3_5.Uxy0_f(j,index_sg(1:3))); anl.y3_5.Uxy0_f(j,index_sg)']*1e-3;
    
    % Uxy0(2:7)=mean(Uxy0(7:10));%
    
    %     %---take mean value
    %            Uxy0=anl.y3_5.Uxy0_f(j,index_sg)*0;
    %            Uxy0=Uxy0 + mean(anl.y3_5.Uxy0_f(j,7:12))*1e-3;
    %            Uxy0=[0;Uxy0'];
    
    
    %G_tmp=[0; G];
    Gamma_Tmp=[0;anl.y3_5.Gamma(index_sg(1)); anl.y3_5.Gamma(index_sg)'];
    
    eqm{j}=Crack_EqMotion(x,Uxy0,Gamma_Tmp);
    Gs_Tmp=interp1(eqm{j}.l*1e3,eqm{j}.Gs,anl.frontRaw{j}.xv);
    Gamma_Tmp=interp1(eqm{j}.l*1e3,eqm{j}.Gamma,anl.frontRaw{j}.xv);
    anl.frontRaw{j}.Gs_Gamma=Gs_Tmp./Gamma_Tmp; %Static energy release rate / Fracture energy
    
    %--- get Cf at sg location
    for l=1:length(index_sg)
        [~, index]=min(abs(anl.frontRaw{j}.xv - anl.y3_5.x_sg(index_sg(l))));
        if (index<4)
        anl.y3_5.Cf(j,index_sg(l))=0;    
        else
        anl.y3_5.Cf(j,index_sg(l))=mean( anl.frontRaw{j}.v(index-3:index+3) );
        end
    end
    %--- get Gs at sg location
    anl.y3_5.Gs(j,:)=interp1(eqm{j}.l*1e3,eqm{j}.Gs,anl.y3_5.x_sg);
    
end

%------------------------------ Plot stuff

% %--------------------plot  K_II^2 and Cf
figN1=figure;
 subplot(1,2,1);
hold all;
 subplot(1,2,2);
hold all;
figN2=figure;
subplot(1,2,1);
hold all;
 subplot(1,2,2);
hold all;
figN3=figure;
hold all;


for j=1:length(anl.e)
    %for j=7:7
    
    figure(figN1);
    subplot(1,2,1);
    plot(eqm{j}.l*1e3,eqm{j}.Gs,'.-');
    subplot(1,2,2);
    plot(eqm{j}.l*1e3,eqm{j}.Uxy0*1e3,'.-');

    
    figure(figN3);
    %plot(eqm{j}.l*1e3,eqm{j}.v.*eqm{j}.Cr,'-');
    plot(eqm{j}.l*1e3,eqm{j}.v,'-');
    
    figure(figN2);
    subplot(1,2,1);
    index1=anl.frontRaw{j}.xv>xStart;
    index2=anl.frontRaw{j}.xv<xEnd;
    index=logical(index1.*index2);
    
    plot(anl.frontRaw{j}.xv(index),smooth(anl.frontRaw{j}.v(index),3)/Cr,'.-')
    subplot(1,2,2);
    plot(anl.frontRaw{j}.Gs_Gamma(index),smooth(anl.frontRaw{j}.v(index),3)/Cr,'.-');
    
end

figure(figN1);
subplot(1,2,1);
legend(num2str(anl.e') );legend off;
title(anl.lgnd{1}(1:end-3));
plot(eqm{j}.l*1e3,eqm{j}.Gamma,'black-');
ylabel('Gs');

subplot(1,2,2);
plot(anl.y3_5.x_sg(index_sg),anl.y3_5.Uxy0_f(:,index_sg),'o');
ylabel('Uxy_0-Uxy_r');
my_legend_add(num2str(anl.e') );legend off;
my_legend_add(num2str(anl.e') );legend off;



figure(figN2);
subplot(1,2,2);
Cf_tmp=linspace(1,Cr,100);
my_legend_add(num2str(anl.e') );legend off
plot(eqm{j}.g(Cf_tmp).^-1,Cf_tmp/Cr,'black-')
title(anl.lgnd{1}(1:end-3));
my_legend_add('g^{-1}' );legend off
xlabel('Gs/G');
ylim([0 eqm{j}.Cr/Cr]);
subplot(1,2,1);
my_legend_add(num2str(anl.e') );legend off
xlabel('x(mm)');
ylim([0 eqm{j}.Cr/Cr]);

figure(figN3);
legend(num2str(anl.e') );legend off
title(anl.lgnd{1}(1:end-3));


% % -----------plot Uii for each event
figN6=figure;
subplot(1,2,1);
hold all;
xlabel('x(mm)');
ylabel('Uxx*Cf(m/s)');
title(anl.lgnd{1}(1:end-3));
subplot(1,2,2);
hold all;
xlabel(' x(mm)');
ylabel('Uyy*Cf(m/s)');

figN7=figure;
subplot(1,2,1);
hold all;
xlabel(' Gs/Gamma');
ylabel('Uxx*Cf(m/s)/Gamma^{0.5}');
ylim([2e-3 2])
title(anl.lgnd{1}(1:end-3));
subplot(1,2,2);
hold all;
xlabel(' Gs/Gamma');
ylabel('Uyy*Cf(m/s)/G^{0.5}');
ylim([2e-3 2]/3)

for j=1:length(anl.e)
    
    %----plot Uii*cf Vs x
    
    figure(figN6);
    subplot(1,2,1);
    plot(anl.y3_5.x_sg(index_sg_Plot),1e-3*anl.y3_5.UxxP(j,index_sg_Plot).*anl.y3_5.Cf(j,index_sg_Plot),'o-');
    my_legend_add( num2str(anl.e(j)) );
    
    subplot(1,2,2);
    plot(anl.y3_5.x_sg(index_sg_Plot),1e-3*anl.y3_5.UyyP(j,index_sg_Plot).*anl.y3_5.Cf(j,index_sg_Plot),'o-');
    my_legend_add( num2str(anl.e(j)) );
    
    % % --------plot Uii*cf Vs K^2
    figure(figN7);
    subplot(1,2,1);
    x=anl.y3_5.Gs(j,index_sg_Plot)./anl.y3_5.Gamma(index_sg_Plot);
    y=1e-3*anl.y3_5.UxxP(j,index_sg_Plot).*anl.y3_5.Cf(j,index_sg_Plot)./anl.y3_5.Gamma(index_sg_Plot).^0.5;
    semilogy(x,y,'o');
    my_legend_add( num2str(anl.e(j)));
    
    subplot(1,2,2);
    y=1e-3*anl.y3_5.UyyP(j,index_sg_Plot).*anl.y3_5.Cf(j,index_sg_Plot)./anl.y3_5.Gamma(index_sg_Plot).^0.5;
    semilogy(x,y,'o');
    my_legend_add( num2str(anl.e(j)));
    
    
end



% %-------------plot Uii Vs Cf  &  Uii*cf Vs K^2 for each sg (all sg)

figN4=figure;
figN5=figure;

for j=1:length(index_sg)
    
    % ---------plot Uii Vs Cf
    
    figure(figN4);
    subplot(1,2,1);
    plot(anl.y3_5.Cf(:,index_sg(j)),anl.y3_5.UxxP(:,index_sg(j)),'o');
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel('Cf');
    ylabel('Uxx');
    title(anl.lgnd{1}(1:end-3));
    xlim([0 Cr]);
    
    subplot(1,2,2);
    plot(anl.y3_5.Cf(:,index_sg(j)),anl.y3_5.UyyP(:,index_sg(j)),'o');
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel('Cf');
    ylabel('Uyy');
    xlim([0 Cr]);
    
    % % --------plot Uii*cf/Gamma^0.5 Vs K^2
    figure(figN5);
    subplot(1,2,1);
    x=anl.y3_5.Gs(:,index_sg(j))./anl.y3_5.Gamma(index_sg(j));
    y=1e-3*anl.y3_5.UxxP(:,index_sg(j)).*anl.y3_5.Cf(:,index_sg(j))./anl.y3_5.Gamma(index_sg(j)).^0.5;
    semilogy(x,y,'o');
    
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel(' Gs/Gamma');
    ylabel('Uxx*Cf(m/s)/Gamma^{0.5}');
    ylim([2e-3 2])
    title(anl.lgnd{1}(1:end-3));
    
    subplot(1,2,2);
    x=anl.y3_5.Gs(:,index_sg(j))./anl.y3_5.Gamma(index_sg(j));
    y=1e-3*anl.y3_5.UyyP(:,index_sg(j)).*anl.y3_5.Cf(:,index_sg(j))./anl.y3_5.Gamma(index_sg(j)).^0.5;
    semilogy(x,y,'o');
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel(' Gs/Gamma');
    ylabel('Uyy*Cf(m/s)/Gamma^{0.5}');
    ylim([2e-3 2]/3)
    
end

%--calculate predicted value.
%load('G:\Frics\Analyze\SuperShear\CohesiveZoneModel.mat');
k=Cs/Cd;%Broberg p.330
SPlus=@(v) CrackSolutionCalcSPlus(1./v,1/Cd,1/Cs);
A=@(v) CrackSolutionCalcA(v,2);
kappa=@(v) 1./SPlus(v).*(1-v/Cr)./(1-v/Cs).^0.5;%Dynamic intensity factor Eq. 6.4.42
g=@(v) kappa(v).^2.*A(v);
GsTmp=linspace(0.05,0.999,20); %%fTmp= @(v) G*4*mu*(1-k^2)-g(v)*KII_s(j)^2;
for j=1:length(GsTmp)
    fTmp=@(v) GsTmp(j)-g(v);
    CfTmp(j)=fzero(fTmp,[0.01 Cr*0.9999])/Cr;
end
%CfTmp=linspace(0.01,0.99,25);
%GTmp=g(CfTmp*Cr);
sol=CrackSolutionForhAnalyze(CfTmp,3.5e-3);

%-------plot predicted values
figure(figN4);
subplot(1,2,1);
plot(sol.v'*sol.Cr,sol.Uxx.max*1e3,'black.-');
my_legend_add(['G=' num2str(sol.Gamma)]);
G_ratio=sol.Gamma/min(anl.y3_5.Gamma);
plot(sol.v'*sol.Cr,sol.Uxx.max*1e3/(G_ratio)^0.5,'black.-');
my_legend_add(['G=' num2str(sol.Gamma/G_ratio)]);

subplot(1,2,2);
plot(sol.v'*sol.Cr,sol.Uyy.max*1e3,'black.-');
my_legend_add(['G=' num2str(sol.Gamma)]);
plot(sol.v'*sol.Cr,sol.Uyy.max*1e3/G_ratio^0.5,'.black-');
my_legend_add(['G=' num2str(sol.Gamma/G_ratio)]);


figure(figN5);
subplot(1,2,1);
%plot((sol.Gamma./GTmp).^1,sol.Uxx.max.*sol.v'*sol.Cr,'.-');
semilogy((1./GsTmp).^1,sol.Uxx.max.*sol.v'*sol.Cr/sol.Gamma^0.5,'black.-');
%plot((sol.Gamma./GTmp).^1,sol.Uxx.max,'.-');
my_legend_add(['G=' num2str(sol.Gamma)]);

subplot(1,2,2);
%plot((sol.Gamma./GTmp).^1,sol.Uyy.max.*sol.v'*sol.Cr,'.-');
semilogy((1./GsTmp).^1,sol.Uyy.max.*sol.v'*sol.Cr/sol.Gamma^0.5,'black.-');
%plot((sol.Gamma./GTmp).^1,sol.Uyy.max,'.-');
my_legend_add(['G=' num2str(sol.Gamma)]);

figure(figN7);
subplot(1,2,1);
semilogy((1./GsTmp).^1,sol.Uxx.max.*sol.v'*sol.Cr/sol.Gamma^0.5,'black.-');
subplot(1,2,2);
semilogy((1./GsTmp).^1,sol.Uyy.max.*sol.v'*sol.Cr/sol.Gamma^0.5,'black.-');
