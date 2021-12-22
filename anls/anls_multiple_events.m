function anl=anls_multiple_events(exp,event,sg_num,slip_smt,baseFig)


%---read and create data
anl.event=event;
for j=1:length(event)
    e=acq132_event_get_data(exp,event(j),'start','end',1,'ch_094_ex','Sxy','Sxx','Syy','x_sg');
    f=calc_Fx_stress_splined(e);
    
    anl.t_trig(j)=e.t_trig;
    anl.t(:,j)=e.t;
    anl.ch2(:,j)=smooth(e.ch_094_ex(:,2),5);
    anl.ch5(:,j)=e.ch_094_ex(:,5);
    anl.ch4(:,j)=e.ch_094_ex(:,4);
    %anl.slip(:,j)=(e.ch_094_ex(:,4)-e.ch_094_ex(:,2));
    anl.slip(:,j)=(e.ch_094_ex(:,5)-e.ch_094_ex(:,2));
    
    anl.velocity=diff(my_smooth(anl.slip,slip_smt))./diff(anl.t)/10;%[cm/sec]
    anl.Sxy(:,j)=e.Sxy(:,sg_num);
    anl.Sxx(:,j)=smooth(e.Sxx(:,sg_num),3);
    anl.Syy(:,j)=smooth(e.Syy(:,sg_num),3);
    
    index=find(f.x==e.x_sg(sg_num))-1;
    anl.dSxxdx(:,j)=f.dSxxdx(:,index);
    anl.Fx(:,j)=f.Fx(:,index);
    
    %---Cut single pix v.s t
    if  exist([exp '\vds'],'dir')==7
        vds=vds_event_get_data(exp,event(j),1,1000,3);
        [~,pix]=min(abs(vds.x-e.x_sg(sg_num)));
        anl.A(:,j)=vds.lines(:,pix);
        anl.A_t(:,j)=vds.t;
        %anl.A=subtruct_norm(vds.lines(:,pix),1); %normilize the line;
    end
    
    if  exist([exp '\Ph'],'dir')==7
        phE=phantomGetLines(exp,event(j),'start','min','end',1000,1);
        A=get_A_at_x(phE,e.x_sg(sg_num),9);
        anl.A(:,j)=A.lines;
        anl.A_t(:,j)=A.t;
    end
    
    
end

if nargin==5  %plot only if baseFig is entered
    
    %---plot
    
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.t,-anl.ch5,'.-');
    % title('up slip');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
    %
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.t,anl.ch3,'.-');
    % title('down slip');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
    
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.t,anl.slip,'.-');
    % title('slip');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
    %
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.t,anl.ch4,'.-');
    % title('corner slip');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
    %
    %
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.t(1:end-1,:)+0.5*diff(anl.t),anl.velocity,'.-');
    % title('velocity');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
    %
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.t,anl.Fx,'.-');
    % title('Fx');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
    %
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.t,anl.Fx./anl.Syy,'.-');
    % title('Fx/Syy');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
    %
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.t,anl.Syy,'.-');
    % title('Syy');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
    %
    % if exist('anl.A','var')==1
    % figure(baseFig);
    % baseFig=baseFig+1;
    % plot(anl.A_t,anl.A,'.-');
    % title('A');
    % legend(cellstr(num2str(anl.event')));
    % legend('off');
end

end
