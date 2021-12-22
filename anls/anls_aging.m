function anl=anls_aging(exper,event)


%-------------
path=pwd;
path=[path(end-9:end) '\' exper '\e'];


phS=phantomGetLines(exper,0,'start','min','end',1,1,'all',2:7);
x_sg=[97.02 105.36 113.02 121.02  128.43];
%x_sg=[105.36 121.02 ];
fig_N=3;

for j=1:length(event)
    
     
   lgnd{j}=[path num2str(event(j))];
    
   %phE=phantomGetLines(exper,event(j),-5e-3,'min',10e-3,1,1,'all',2:7);
   %line1=mean(phE.lines(1:100,:),1);
   %phE.lines=phE.lines./repmat(line1,length(phE.t),1);
   
   for k=1:length(x_sg)
       
       
       [~,index_x]=min(abs(phS.x-x_sg(k)));
       dx=20;
       
      % index_t=find(subtruct_norm(phE.lines(:,index_x),1)<0.95,1,'first');
      % t_event=phE.t(index_t);
       
       figure(fig_N-1+k);
       
       %--fast
      % tf=phE.t-t_event;
      % Af=mean(phE.lines(:,index_x-dx:index_x+dx),2);
      % Af=Af(tf>0);
      % tf=tf(tf>0);
       
       %----slow 
       
%       ts=phS.t-phE.t_trigger-t_event;
       ts=phS.t-t_event;

       As=mean(phS.lines(:,index_x-dx:index_x+dx),2);
       As=As./mean(line1(index_x-dx:index_x+dx));
       As=As(ts>0);
       ts=ts(ts>0);
       
       %---combine slow and fast
       A=[ As];
       t=[ ts];
       %A=[Af; As];
       %t=[tf; ts];
       
             
              
       semilogx(t,A,'.-');
       hold all;
       title([phS.Date '/' num2str(phS.Exp) 'x= ' num2str(x_sg(k))]);
       my_legend_add(num2str(event(j)));
   end
   


    
    
end

% anl.lgnd=lgnd;
% anl.event=event;
% anl.event_slow=e_Slow;
% anl.front=front;
% anl.phECut=phECut;