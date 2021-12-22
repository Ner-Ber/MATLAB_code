function anl=anls_aging2(exper,eventV)

%Only from slow measurements of A

%-------------
path=pwd;
path=[path(end-9:end) '\' exper '\e'];


phS=phantomGetLines(exper,0,'start','min','end',1,1,'all',2:7);
%x_sg=[89.45 97.02 105.36 113.02 121.02  128.43 135.9];
x_sg=[105.36 113.02 121.02  128.43 135.9];

dx=20;

fig_N=figure;
fig_N=fig_N.Number;

%------Find t0
[~,index_x]=min(abs(phS.x-x_sg(1)));
As=mean(phS.lines(:,index_x-dx:index_x+dx),2);
D_As=diff(As);
index=1:length(phS.t);
index_t0=index(-D_As>150);
index_tend=index_t0(2:end);
t0=(phS.t(index_t0)+phS.t(index_t0+1))/2;

for j=1:length(eventV)
    
    event=eventV(j);
    lgnd{j}=[path num2str(event(j))];
    
    figure(fig_N+j-1);
    
    for k=1:length(x_sg)
        
        
        [~,index_x]=min(abs(phS.x-x_sg(k)));
        
        
        A=mean(phS.lines(:,index_x-dx:index_x+dx),2);
        A0=mean(A(index_t0-10:index_t0-8));
        t=phS.t(index_t0(event):index_tend(event))-t0(event(j));
        A=A(index_t0(event):index_tend(event));
        A=A/A0;
        
        %--cut after 100s
        A=A(t<100);
        t=t(t<100);
        t=t(3:end);
        A=A(3:end);
        semilogx(t(1:end),A(1:end),'.-');
%        hold all;
        %my_legend_add(num2str(x_sg(k)));
        
        %-fit
        index=round(logspace(0,log10(length(t))));
        A=A(index);
        t=t(index);
        %semilogx(t(1:end),A(1:end),'.-');
                hold all;
        t=log(t);
        
        f=polyfit(t,A,1);
        my_legend_add([num2str(x_sg(k)),'b=' num2str(f(1)) 'c=' num2str(f(2))]);
        x=logspace(-2,3,100);
        plot(x,f(2)+f(1)*log(x),'k');
        my_legend_add([num2str(x_sg(k)),'b=' num2str(f(1)) 'c=' num2str(f(2))]);

    end
    
    title([phS.Date '/' num2str(phS.Exp) 'e= ' num2str(event) 't=' num2str(t0(event(j)))]);
    
    
    
    
end

% anl.lgnd=lgnd;
% anl.event=event;
% anl.event_slow=e_Slow;
% anl.front=front;
% anl.phECut=phECut;