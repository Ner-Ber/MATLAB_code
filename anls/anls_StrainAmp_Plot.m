function anls_StrainAmp_Plot(anl,kind)

%gets analEqm and plots the amplitudes. 

[Cd Cs Cr nu ro E mu Gamma PlaneStrain tau_p Xc0]=CrackSolutionMaterialProperties;


%index_sg=7:(length(anl.y3_5.x_sg)); %Choose sg for plot
index_sg=8:14; %Choose sg for plot

%index_e=[1:11 13:length(anl.e)];
index_e=1:length(anl.e);

if (length(index_e)<length(anl.e))
    display('Not all event are plotted')
end
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


%-------------plot Uii Vs Cf for each sg

figN1=figure;
figN2=figure;
figN3=figure;

for j=1:length(index_sg)
    
        
    % ---------plot Uii Vs Cf
    
    %---prepare color based on Uxy0
    %c=anl.y3_5.Uxy0_f*1e-3*2*mu/tau_p;
    
    
    figure(figN1);
    plot(anl.y3_5.Cf(index_e,index_sg(j)),anl.y3_5.UxxP(index_e,index_sg(j)),'o');
    %scatter(anl.y3_5.Cf(index_e,index_sg(j)),anl.y3_5.UxxP(index_e,index_sg(j)),[],c(index_e,index_sg(j)));
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel('Cf');
    ylabel('Uxx');
    title(anl.lgnd{1}(1:end-3));
    %xlim([0 Cr*1.1]);
    
    figure(figN2);
    plot(anl.y3_5.Cf(index_e,index_sg(j)),anl.y3_5.UyyP(index_e,index_sg(j)),'o');
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel('Cf');
    ylabel('Uyy');
    %xlim([0 Cr*1.1]);
    
    figure(figN3);
    plot(anl.y3_5.Cf(index_e,index_sg(j)),anl.y3_5.Uxy0_f(index_e,index_sg(j)),'o');
    hold all;
    my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    xlabel('Cf');
    ylabel('Uxy');
    
    
%     if strcmp(kind,'Rayleigh')
%         %Sub-Rayleigh
%         tau_r=calc_dynamic_stress(anl.y3_5.Uxyf3(index_e,index_sg(j)),anl.y3_5.Uxy0(index_e,index_sg(j)));
%         
%     else
%         %SuperShear
%         tau_r=calc_dynamic_stress(anl.y3_5.Uxyf(index_e,index_sg(j)),anl.y3_5.Uxy0(index_e,index_sg(j)));
%         
%     end
  
end


% % %-------------plot Uii Vs Cf for each event
% 
% figN1=figure;
% figN2=figure;
% %figN3=figure;
% 
% for j=1:length(index_e)
%     
%         
%     % ---------plot Uii Vs Cf
%     
%     figure(figN1);
%     plot(anl.y3_5.Cf(index_e(j),index_sg),anl.y3_5.UxxP(index_e(j),index_sg),'o');
%     hold all;
%     my_legend_add( num2str(anl.e(index_e(j))) );
%     xlabel('Cf');
%     ylabel('Uxx');
%     title([anl.lgnd{1}(1:end-3) '--x=' num2str(anl.y3_5.x_sg(index_sg))]);
%     %xlim([0 Cr*1.1]);
%     
%     figure(figN2);
%     plot(anl.y3_5.Cf(index_e(j),index_sg),anl.y3_5.UyyP(index_e(j),index_sg),'o');
%     hold all;
%     my_legend_add( num2str(anl.e(index_e(j))) );
%     xlabel('Cf');
%     ylabel('Uyy');
%     %xlim([0 Cr*1.1]);
%     
%     
%     
% %     if strcmp(kind,'Rayleigh')
% %         %Sub-Rayleigh
% %         tau_r=calc_dynamic_stress(anl.y3_5.Uxyf3(index_e,index_sg(j)),anl.y3_5.Uxy0(index_e,index_sg(j)));
% %         
% %     else
% %         %SuperShear
% %         tau_r=calc_dynamic_stress(anl.y3_5.Uxyf(index_e,index_sg(j)),anl.y3_5.Uxy0(index_e,index_sg(j)));
% %         
% %     end
%   
% end
