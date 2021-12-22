function anls_Af_Plot(anl,kind,solSR)

%After using anls_A3 and anls_EquationOfMotion plot Ar Vs Uxx*Cf

index_sg=1:(length(anl.y3_5.x_sg)); %Choose sg for plot
%index_sg=11;
%index_e=[1 2:length(anl.event)];
index_e=[1 2 3 4:length(anl.event)];
anl.y3_5.Cf=zeros(length(anl.event),length(anl.y3_5.x_sg));

for j=1:length(anl.event)
    
       
        if strcmp(kind,'Rayleigh')
        %SubRayleigh
        %anl.frontRaw{j}.v(anl.frontRaw{j}.v>1300)=NaN;
        anl.frontRaw{j}.v(anl.frontRaw{j}.v>1300)=NaN;
    else
        %SuperShear
        anl.frontRaw{j}.v(anl.frontRaw{j}.v<1300)=NaN;
        end
    
        %--- get Cf at sg location
     
    for l=1:length(index_sg)
        
      [~, index1]=min(abs(anl.frontRaw{j}.xv - (anl.y3_5.x_sg(index_sg(l))-12)));
        [~, index2]=min(abs(anl.frontRaw{j}.xv - (anl.y3_5.x_sg(index_sg(l))+12)));
        
        anl.y3_5.Cf(j,index_sg(l))=mean( anl.frontRaw{j}.v(index1:index2) );
    
        %--- get Af at sg location
        xStart=anl.y3_5.x_sg(index_sg(l))-5;
        xEnd=anl.y3_5.x_sg(index_sg(l))+5;
        %---For each strain gage take average over +-5mm
        
        %anl.y3_5.Af(j,index_sg(l))=mean(anl.front{j}.Af((anl.front{j}.x1>xStart&anl.front{j}.x1<xEnd)));
        %anl.y3_5.Af_from_x(j,index_sg(l))=mean(anl.phECut{j}.Af_from_x(anl.phECut{j}.frontX>xStart & anl.phECut{j}.frontX<xEnd) );
        %anl.y3_5.Af_from_t(j,index_sg(l))=mean(anl.phECut{j}.Af_from_t(anl.phECut{j}.x>xStart & anl.phECut{j}.x<xEnd) );
        
        anl.y3_5.Af_from_t(j,index_sg(l))=mean(anl.phECut{j}.Af3_from_t(anl.phECut{j}.x>xStart & anl.phECut{j}.x<xEnd) );
        %A0(j,index_sg(l))=mean(anl.phECut{j}.firstline((anl.phECut{j}.x>xStart & anl.phECut{j}.x<xEnd)));  
    end
    
    A0(j)=mean(anl.phECut{j}.firstline((anl.phECut{j}.x>90 & anl.phECut{j}.x<140)));%100-130

end



% %-------------plot Af Vs Uxx*Cf for each sg

%figN1=figure;
%figN2=figure;
figN3=figure;

for j=1:length(index_sg)
    
    % ---------plot Ar Vs -Uxx*Cf
    
    Cf=anl.y3_5.Cf(index_e,index_sg(j));
    
    %UxxP=anl.y3_5.UxxP(index_e,index_sg(j))*1e-3;
    %v=UxxP.*Cf;
    %v=v*2*0.5;
    %v=calc_vSlip(v,solSR);
    
    Uxxf=anl.y3_5.Uxxf3(index_e,index_sg(j))*1e-3;
    v=2*Uxxf.*Cf;
       
     v(v>0.15)=NaN;
    %plot(UxxP.*Cf,Uxxf.*Cf,'o');
    %hold all;
    
%         figure(figN1);
%         Af=anl.y3_5.Af(index_e,index_sg(j));
%         semilogx(v,Af,'o');
%         hold all;
%         my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    
    %     figure(figN2);
    %     Af=anl.y3_5.Af_from_x(index_e,index_sg(j));
    %     semilogx(v,Af,'o');
    %     hold all;
    %     my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))) );
    
      %figure(figN3);
      Af=anl.y3_5.Af_from_t(index_e,index_sg(j));
         
  %Af=(anl.y3_5.Af_from_t(index_e,index_sg(j)).*A0(index_e)'-320)/(1280-320);
  %Af=(anl.y3_5.Af_from_t(index_e,index_sg(j)).*A0(1,index_sg(j))-320)/(1290-320);
  
   %Af=(anl.y3_5.Af_from_t(index_e,index_sg(j))*1290-320)/(1290-320);
       %Af=Af*0.98.*(1+0.03*log(1+v/0.3));
       
      semilogx(v,Af,'k.');
      hold all;
      my_legend_add( num2str(anl.y3_5.x_sg(index_sg(j))));
    
end
% figure(figN1);
% %xlabel('Uxx*Cf(m/s)');
% ylabel('Af/A0');
% xlim([5e-3 2]);
% ylim([0.6 0.9]);
% title([anl.lgnd{1}(1:end-3) ' Af_from_x']);

%figure(figN3);
%xlabel('Uxx*Cf(m/s)');
%ylabel('Af/A0');
%xlim([5e-3 2]);
%ylim([0.6 0.9]);
%title([anl.lgnd{1}(1:end-3) ' Af_from_t'] );

% 
% figure;
% Gamma=1.2;
%  v3_5=solSR2{2}.Uxx.max.*solSR2{2}.v'*solSR2{2}.Cr*(Gamma/solSR2{1}.Gamma)^0.5;
%  v_0=solSR2{1}.Uxx.max.*solSR2{1}.v'*solSR2{1}.Cr*(Gamma/solSR2{1}.Gamma)^0.5;
% 
% 
% for j=1:length(fig.x) vSlip{j}=interp1(v3_5,2*v_0,fig.x{j}); end
% for j=1:length(fig.x) semilogx(vSlip{j},fig.y{j},'ko'); hold all; end
% legend(fig.DisplayName);legend off;
% %ylim([0.5 0.9]);
% xlim([1e-2 30]);
% a=get(gca,'Children');
% for j=1:length(a) set(a(j),'color', fig.color{length(a)-j+1}); end   

