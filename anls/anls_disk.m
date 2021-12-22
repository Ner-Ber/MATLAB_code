function acqSlop=anls_disk(exper)


%for j=1:length(exper) acq132_stream_save_mat(exper{j},'stream'); end

for j=1:length(exper)
    acqSTmp=acq132_stream_load_mat([exper{j} '\stream.mat'],'start','min','end',10,'N');
    [~,j2]=max(acqSTmp.N);
    [~,j1]=min(abs(acqSTmp.N(1:j2)-50));
    startT=acqSTmp.t(j1);
    endT=acqSTmp.t(j2);
    clear acqSTmp
    acqS{j}=acq132_stream_load_mat([exper{j} '\stream.mat'],startT,'min',endT,10,'N','U1','U2','U3','sg_angle');
    
    %[acqS{j}.U1,acqS{j}.U2,acqS{j}.U3]=calc_shear_sensitivity4(acqS{j}.U1,acqS{j}.U2,acqS{j}.U3,[0.07,0.1,0.9,-0.04]);
    
    %[acqS{j}.U1,acqS{j}.U2,acqS{j}.U3]=calc_shear_sensitivity4(acqS{j}.U1,acqS{j}.U2,acqS{j}.U3,[0.027,0.1,0.91,-0.08]);
    %[acqS{j}.U1,acqS{j}.U2,acqS{j}.U3]=calc_shear_sensitivity4(acqS{j}.U1,acqS{j}.U2,acqS{j}.U3,[0.079,0.1,0.9,-0.03]);
    
    [acqS{j}.U1,acqS{j}.U2,acqS{j}.U3]=calc_shear_sensitivity4(acqS{j}.U1,acqS{j}.U2,acqS{j}.U3,[0,00,0.95,-0.08]);
    
    [~,~,~,acqS{j}.Uxx,acqS{j}.Uyy,acqS{j}.Uxy,~]=calculate_stress_strain(acqS{j}.U1,acqS{j}.U2,acqS{j}.U3,acqS{j}.sg_angle);
    
    %linear fit
    for k=1:2
        
        a=fit(acqS{j}.N,acqS{j}.U1(:,k),'poly1');
        acqSlop{k}.U1(1,j)=a.p1;
        
        a=fit(acqS{j}.N,acqS{j}.U2(:,k),'poly1');
        acqSlop{k}.U2(1,j)=a.p1;
        
        a=fit(acqS{j}.N,acqS{j}.U3(:,k),'poly1');
        acqSlop{k}.U3(1,j)=a.p1;
        
        acqSlop{k}.sg_angle(j)=acqS{j}.sg_angle(k);
    end
    
end

%---Theory

x=(50:0.01:150);
Uxx=2*x'/3e6/pi/7.6e-3/50e-3*9.8;
Uyy=-3.3*x'/3e6/pi/7.6e-3/50e-3*9.8;
Uxx=Uxx/1.08;
Uyy=Uyy/1.08;
[Uxx,Uxy,Uyy]=Rotate_strain(Uxx,Uxx*0,Uyy,-2/180*pi);
[~,~,U45]=Rotate_strain(Uxx,Uxy,Uyy,45/180*pi);
[~,~,U_45]=Rotate_strain(Uxx,Uxy,Uyy,-45/180*pi);


% %------Plot
% figure;
% hold all;
% k=2;
% for j=1:length(exper)
% plot(acqS{j}.N,[subtruct_norm(acqS{j}.U1(:,k)) subtruct_norm(acqS{j}.U2(:,k) ) subtruct_norm( acqS{j}.U3(:,k)) subtruct_norm( acqS{j}.U1(:,k)+acqS{j}.U3(:,k)-acqS{j}.U2(:,k) ) ]);
% title(['side' num2str(k)]);
% xlim([50 150]);
% ylim([-1.25 0.8]);
%
% end
% plot(x,[subtruct_norm(Uyy) subtruct_norm(Uxx) subtruct_norm(U45) subtruct_norm(U_45)],'black');


% figure;
% hold all;
% k=2;
% for j=1:length(exper)
% plot(acqS{j}.N,subtruct_norm([acqS{j}.Uxx(:,k) acqS{j}.Uxy(:,k) -acqS{j}.Uyy(:,k)]));
% title(['side' num2str(k)]);
% xlim([50 150]);
% ylim([-1.25 0.8]);
% end
% plot(x,[subtruct_norm(Uxx) subtruct_norm(Uxy) -subtruct_norm(Uyy)],'black');

figure;
hold all;
for k=1:2
for j=1:length(exper)
    plot(acqS{j}.N,subtruct_norm(acqS{j}.U1(:,k)) + subtruct_norm( acqS{j}.U3(:,k)) );
    title(['side' num2str(k)]);
    xlim([50 150]);
    ylim([-1.25 0.8]);
end
end
plot(x,subtruct_norm(Uyy)+subtruct_norm(Uxx) ,'black');



% figure;
% hold all;
% for j=1:length(exper)
% plot(acqS{j}.N,subtruct_norm(acqS{j}.Uxx(:,1) + acqS{j}.Uyy(:,1)));
% ylabel('Uxx+Uyy');
% title('side1');
% my_legend_add(num2str(acqS{j}.sg_angle(1)/pi*180) );
% xlim([50 150]);
% ylim([-0.6 0.1]);
% end
%
% figure;
% hold all;
% for j=1:length(exper)
% plot(acqS{j}.N,subtruct_norm([acqS{j}.Uxx(:,1) acqS{j}.Uyy(:,1) acqS{j}.Uxy(:,1)]));
% title('side1');
% xlim([40 150]);
% ylim([-1.25 0.8]);
% end
%
% figure;
% hold all;
% for j=1:length(exper)
% plot(acqS{j}.N,subtruct_norm([acqS{j}.Uxx(:,2) + acqS{j}.Uyy(:,2)]));
% ylabel('Uxx+Uyy');
% title('side2');
% my_legend_add(num2str(acqS{j}.sg_angle(2)/pi*180) );
% xlim([50 150]);
% ylim([-0.6 0.1]);
% end
%
% figure;
% hold all;
% for j=1:length(exper)
% plot(acqS{j}.N,subtruct_norm([acqS{j}.Uxx(:,2) acqS{j}.Uyy(:,2) acqS{j}.Uxy(:,2)]));
% title('side2');
% xlim([40 150]);
% ylim([-1.25 0.8]);
% end
