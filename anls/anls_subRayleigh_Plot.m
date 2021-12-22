function anls_subRayleigh_Plot(anl)

[Cd Cs Cr, ~ , ~ ,E mu Gamma, ~ ,~]=CrackSolutionMaterialProperties;



index_sg=2:(length(anl.y3_5.x_sg)); %Choose sg for plot
index_e=[1:4 6:length(anl.e)];
anl.y3_5.Cf=zeros(length(anl.e),length(anl.y3_5.x_sg));

for j=1:length(anl.e)
    %--- get Cf at sg location
    
    for l=1:length(index_sg)
        
        [~, index]=min(abs(anl.frontRaw{j}.xv - anl.y3_5.x_sg(index_sg(l))));
        if (index<4)
            index=4;
        end
        anl.y3_5.Cf(j,index_sg(l))=mean( anl.frontRaw{j}.v(index-3:index+3) );
    end
end


% %-------------plot Uii Vs Cf for each sg 

figN1=figure;

for j=1:length(index_sg)
    
    % ---------plot Uii Vs Cf
    
    figure(figN1);
    subplot(1,2,1);
    plot(anl.y3_5.Cf(index_e,index_sg(j)),anl.y3_5.UxxP(index_e,index_sg(j)),'o');
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel('Cf');
    ylabel('Uxx');
    title(anl.lgnd{1}(1:end-3));
    xlim([0 Cr]);
    ylim([0 0.6]);
    
    subplot(1,2,2);
    plot(anl.y3_5.Cf(index_e,index_sg(j)),anl.y3_5.UyyP(index_e,index_sg(j)),'o');
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel('Cf');
    ylabel('Uyy');
    xlim([0 Cr]);
    ylim([0 0.6]/3);
end

% %-------------plot Uii Vs Cf for event 

figN2=figure;

for j=1:length(index_e)
    
    % ---------plot Uii Vs Cf
    
    figure(figN2);
    subplot(1,2,1);
    plot(anl.y3_5.Cf(index_e(j),index_sg),anl.y3_5.UxxP(index_e(j),index_sg),'o');
    hold all;
    my_legend_add( num2str(anl.e(index_e(j))) );
    xlabel('Cf');
    ylabel('Uxx');
    title(anl.lgnd{1}(1:end-3));
    xlim([0 Cr]);
    ylim([0 0.6]);
    
    subplot(1,2,2);
    plot(anl.y3_5.Cf(index_e(j),index_sg),anl.y3_5.UyyP(index_e(j),index_sg),'o');
    hold all;
    my_legend_add( num2str(anl.e(index_e(j))) );
    xlabel('Cf');
    ylabel('Uyy');
    xlim([0 Cr]);
    ylim([0 0.6]/3);
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


figure(figN1);
subplot(1,2,1);
plot(sol.v'*sol.Cr,sol.Uxx.max*1e3,'black.-');
my_legend_add(['G=' num2str(sol.Gamma)]);
G_ratio=sol.Gamma/1.67;
plot(sol.v'*sol.Cr,sol.Uxx.max*1e3/(G_ratio)^0.5,'black.-');
my_legend_add(['G=' num2str(sol.Gamma/G_ratio)]);

subplot(1,2,2);
plot(sol.v'*sol.Cr,sol.Uyy.max*1e3,'black.-');
my_legend_add(['G=' num2str(sol.Gamma)]);
plot(sol.v'*sol.Cr,sol.Uyy.max*1e3/G_ratio^0.5,'.black-');
my_legend_add(['G=' num2str(sol.Gamma/G_ratio)]);

figure(figN2);
subplot(1,2,1);
plot(sol.v'*sol.Cr,sol.Uxx.max*1e3,'black.-');
my_legend_add(['G=' num2str(sol.Gamma)]);
G_ratio=sol.Gamma/1.67;
plot(sol.v'*sol.Cr,sol.Uxx.max*1e3/(G_ratio)^0.5,'black.-');
my_legend_add(['G=' num2str(sol.Gamma/G_ratio)]);

subplot(1,2,2);
plot(sol.v'*sol.Cr,sol.Uyy.max*1e3,'black.-');
my_legend_add(['G=' num2str(sol.Gamma)]);
plot(sol.v'*sol.Cr,sol.Uyy.max*1e3/G_ratio^0.5,'.black-');
my_legend_add(['G=' num2str(sol.Gamma/G_ratio)]);