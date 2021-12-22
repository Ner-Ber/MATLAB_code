function anls_taur_Plot(anl,kind,solSR)

[Cd Cs Cr, ~ , ~ ,E mu Gamma, ~ ,~]=CrackSolutionMaterialProperties;


index_sg=3:(length(anl.y3_5.x_sg)); %Choose sg for plot
index_e=[1 3:length(anl.e)];
%index_e=[12 13];

anl.y3_5.Cf=zeros(length(anl.e),length(index_sg));

for j=1:length(anl.e)
    
    if strcmp(kind,'Rayleigh')
        %SubRayleigh
        anl.frontRaw{j}.v(anl.frontRaw{j}.v>1300)=NaN;
    else
        %SuperShear
        anl.frontRaw{j}.v(anl.frontRaw{j}.v<1300)=NaN;
    end
    
    
    %--- get Cf at sg location
    
    for l=1:length(index_sg)
        
        [~, index1]=min(abs(anl.frontRaw{j}.xv - (anl.y3_5.x_sg(index_sg(l))-12)));%12
        [~, index2]=min(abs(anl.frontRaw{j}.xv - (anl.y3_5.x_sg(index_sg(l))+12)));
        
        anl.y3_5.Cf(j,index_sg(l))=mean( anl.frontRaw{j}.v(index1:index2) );
    end
end


% ------ Compute Uxyf and Uxy0
%
% for j=1:length(anl.y3_5.x_sg)
%     [~,index]=min(abs(anl.pre.x_sg-anl.y3_5.x_sg(j)));
%
%     anl.y3_5.Uxy0(:,j)=anl.pre.Uxy(:,index);
% end

anl.y3_5.Uxyf=anl.y3_5.Uxy0 - anl.y3_5.Uxy0_f;

%--calc Syy
strains.Uyy=anl.y3_5.Uyy0;
strains.Uxx=anl.y3_5.Uxx0;
strains.Uxy=anl.y3_5.Uxy0;
[Syy,~,~]=calcStressStaticFromStrain(strains,'false');
anl.y3_5.Syy=mean(Syy,1);
clear Syy
% %-------------plot Uii Vs Cf for each sg

%figN1=figure;
figN2=figure;
%figN3=figure;
%tau_r0=[1.70,1.94,1.84,2.21,2.45,2.41,2.42,2.42,2.52,2.48,2.64,2.48,2.37,2.32];%2016-05-27\12-12-30

for j=1:length(index_sg)
    
    Cf=anl.y3_5.Cf(index_e,index_sg(j));
    
    if strcmp(kind,'Rayleigh')
        %Sub-Rayleigh
        tau_r=calc_dynamic_stress(anl.y3_5.Uxyf3(index_e,index_sg(j)),anl.y3_5.Uxy0(index_e,index_sg(j)));
        
    else
        %SuperShear
        tau_r=calc_dynamic_stress(anl.y3_5.Uxyf(index_e,index_sg(j)),anl.y3_5.Uxy0(index_e,index_sg(j)));
        
    end
    
    %-------velocity calculation
    %     UxxP=anl.y3_5.UxxP(index_e,index_sg(j))*1e-3;
    %     v=UxxP.*Cf;
    %     v=calc_vSlip(v,solSR);
    
    Uxxf=anl.y3_5.Uxxf3(index_e,index_sg(j))*1e-3;
    v=2*Uxxf.*Cf;
    
    %v(v>0.15)=NaN;
         
    tau_r0(j)=mean(tau_r(v<0.105));
    %tau_r0=max(tau_r);
    
    %     %--calc Syy
    %     strains.Uyy=anl.y3_5.Uyy0(index_e,index_sg(j));
    %     strains.Uxx=anl.y3_5.Uxx0(index_e,index_sg(j))*0;
    %     strains.Uxy=anl.y3_5.Uxy0(index_e,index_sg(j));
    %     [Syy,~,~]=calcStressStaticFromStrain(strains,'false');
    
    Syy=anl.y3_5.Syy(index_sg(j));
    % ---------plot tau_r Vs -Uxx*Cf
%         figure(figN1);
%            semilogx(v,tau_r,'o');
%            hold all;
%            my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
% %    xlim([1e-2 30]);
%         xlabel('v(m/s)');
    
     figure(figN2);
%     semilogx(v,tau_r./tau_r0(j),'o');
%     my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
%     hold all;
%     %xlim([1e-2 30]);
%     xlabel('v(m/s)');
    
    semilogx(v,tau_r./Syy,'o');
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    hold all;
    %     xlim([5e-3 2]);
    
    %     plot(tau_r*0+Syy,tau_r,'o');
    %     my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    %     hold all;
    
    %     figure(figN3);
    %     plot(Cf,UxxP,'o');
    %     my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    %     hold all;
    %
  
end

% figure(figN1);
% xlabel('Uxx*Cf(m/s)');
% ylabel('tau_r(Mpa)');
% title(anl.lgnd{1}(1:end-3));

% figure(figN2);
% ylabel('tau_r/Syy');
% title(anl.lgnd{1}(1:end-3));
