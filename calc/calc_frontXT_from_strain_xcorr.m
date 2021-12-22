function t_front=calc_frontXT_from_strain_xcorr(struct)


smt=1;
struct.Uxx=my_smooth(subtruct_norm(struct.Uxx,2),smt);
Uxx=-(struct.Uxx(:,1:end));
x_sg=struct.x_sg;
t=struct.t;
dt=mean(diff(t));
t_front=0;
t_front_index=1;

%find shift and correction betwean sequent sg.
for j=10:length(x_sg)-1
    
    [~,index]=max(Uxx(:,j));
    t_max1=t(index);
    indexStart1=index-20;
    indexEnd1=index+25;
    
    [~,index]=max(Uxx(:,j+1));
    t_max2=t(index);
    indexStart2=index-20;
    indexEnd2=index+25;
        
    [c,lags]=xcorr(Uxx(indexStart1:indexEnd1,j),Uxx(indexStart2:indexEnd2,j+1),'coeff');%It's better to cut before xcorr
    [~,index]=max(c);
    shift=lags(index);
    t_front(t_front_index+1)=t_front(t_front_index)+t_max2-t_max1-dt*shift;
    t_front_index=t_front_index+1;
    
    %Uxx_shifted=(circshift(Uxx(:,j+1),shift(j)));
%     figure(10);plot(struct.t-t,Uxx(:,j),'.-');hold all
%     figure(10);plot(struct.t-t_front(j+1),[Uxx(:,j+1),Uxx_shifted],'.-');
%      xlim([-0.05 0.05])
%      hold off
    %figure(10);plot(struct.t,Uxx_shifted1,'.-');hold all;
end




