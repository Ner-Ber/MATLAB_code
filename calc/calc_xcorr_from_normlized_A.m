function tlag=calc_xcorr_from_normlized_A(t,A,tStart,tEnd)
%the function is useful if you want to cross corelate functions that are
%different in length. for example if A is in different locations and with
%different front velocities.

%t- matrix of time 
%A- matrix of normlized contact area. 
%tStart,tEnd the time window for crosscor



dt=mean(diff(t),1);
dt=min(dt);
tt=tStart:dt:tEnd;
p=0.9;

tlag(1)=0;
for j=2:length(A(1,:))
    
%correlation is found fron diff of the interplation
    aa1=csaps(t(:,j-1),A(:,j-1),p,tt);
    aa1=diff(aa1);
    aa2=csaps(t(:,j),A(:,j),p,tt);
    aa2=diff(aa2);
    
    [C,lag]=xcorr(aa1,aa2,'biased');
    [~,index]=max(C);
    lag_max_Corr=lag(index);
    tlag(j)=dt*lag_max_Corr;
    
%     %----plot for check
%     figure(3); 
%     plot(dt*lag,C,'.-')
%     
     
%     figure(4);
%    % plot(t(:,j-1),A(:,j-1),'.-');
%      hold all;
%     % plot(tt,aa1,'.-');
%      %plot(t(:,j),A(:,j),'.-');
%      plot(t(:,j)+tlag(j),A(:,j),'.-');
%      %plot(tt+tlag,aa2,'.-');
%      xlim([2*tStart 2*tEnd]);
%     % hold off

end

tlag=cumsum(tlag);
