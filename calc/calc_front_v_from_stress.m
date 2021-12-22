function [sg t_front pk_front v_sg]=calc_front_v_from_stress(t_acq,x_sg,Sxy)

tresh_max=0.016;
tresh_min=0.013;

%[~,t_acq,x_sg,~,~,~,~,~,Sxy,~,~]=acq132_event_get_data(exp_dir,event,pre,post,smt,{'Sxy','Sxx','Syy'});

FirstLine=repmat(mean(Sxy(1:20,:)),length(t_acq),1);
D=Sxy./FirstLine-1;
sg=1:15;

%-----check sg had a drop
for sg_num=1:15
   % if min(D(:,sg_num))<tresh_min
    if mean(D(1:20,sg_num))<mean(D(end-20:end,sg_num))
        sg(sg_num)=0;
    end
end

%----find maximun that exceeds threshold
sg=nonzeros(sg); %sg's that had stress drop
[pk_stress,pk_index]=max(D(:,sg));
[~,sg_tresh]=find(pk_stress>tresh_max); %find those sg that over treshold,sg_tresh is not sg vector in general
sg_no_pk=sg;
sg_no_pk(sg_tresh)=0;
sg_no_pk=nonzeros(sg_no_pk);

x_front=x_sg(sg(sg_tresh));
t_front=t_acq(pk_index(sg_tresh));
pk_front=pk_stress(sg_tresh);

%--------take care of those without peak (sg_no_pk)

for j=1:length(sg_no_pk)
t_index=find(D(:,sg_no_pk(j))<-tresh_min,1);
%[~,t_index]=min(diff(D(:,sg_no_pk(j))));
t_frontt(j,1)=t_acq(t_index);
pk_frontt(j)=D(t_index,sg_no_pk(j));
end

x_frontt=x_sg(sg_no_pk);
[x_front,sort_index]=sort([x_frontt x_front]);
t_front=[t_frontt ; t_front]';
t_front=t_front(sort_index);
pk_front=[pk_frontt pk_front];
pk_front=pk_front(sort_index);

%---calculate front velocity
v=diff(x_front)./diff(t_front);
x_v=x_front(1:end-1)+diff(x_front)/2;
v_sg=spline(x_v,v,x_front); %v at sg position


figure;
plot(x_v,v,'.');
hold all
plot(x_front,v_sg,'.');
figure;
plot(t_acq,D,'.-',t_front,pk_front,'blacko');



