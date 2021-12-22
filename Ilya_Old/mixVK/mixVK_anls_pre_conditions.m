function pre=mixVK_anls_pre_conditions(exp,event,Tstart,lineNum,smtXsg)
%event - Vector
%Tstart - [msec] . mean vaule is taken over +-0.5 msec
%The function calculates all pre event initial conditions.
% 'F','N','Uxy','Uxx','Uyy','Sxy','Sxx','Syy','lines','A','Fx'
% The data is struct that contains matrix for each variable
    
for k=1:length(event)
    
    e=mixVK_acq132_event_get_data(exp,event(k),Tstart-0.5,Tstart+0.5,1,'x_sg','F','N','Uxy','Uxx','Uyy','Sxy','Sxx','Syy');
    %e=smooth_sg_x(e,smtXsg); % could be done smarter and faster
    phE=phantomGetLines(exp,event(k),Tstart-0.5,'min',Tstart+0.5,1000,1,'',lineNum);
    e.lines=phE.lines;
    e.x=phE.x;

    if(k==1)
        %find fields that change
        fNames=fieldnames(e);
        index=find(ismember(fNames, [{'x_sg'},{'exp'},{'x'},{'Date'},{'t'},{'note'}])==0);
        fNames=fNames(index);
        
        %take mean vaule
        for j=1: length(fNames) %the first one is x_sg
            e.(fNames{j})=mean(e.(fNames{j}));
        end
        pre=e;
    else
        %take mean vaule and join the data to previous events
        for j=1: length(fNames)
            pre.(fNames{j})=[pre.(fNames{j}); mean(e.(fNames{j}))];
        end
    end
    
end

%pre.Sxx(6)=mean(pre.Sxx([5,7]));
%pre.Sxx(6)=mean(pre.Sxy([5,7]));

pre.Fx=calc_Fx_stress_splined(pre);
pre.A=get_A_at_x(pre,pre.x_sg,3);
pre.mu=pre.Sxy./pre.Syy;
pre.muA=pre.Sxy./pre.A.lines;
pre.x_ph=pre.x;

pre=rmfield(pre,'x');
