function f=calc_frontXT_from_CursorsSg(x_sg)
%by using cursors on fig eith time seris of sg.
%1.put cursor on the fig
%2.run the function with sg position vector


%---get the cursors 
dcm_obj = datacursormode(gcf);
cursors = getCursorInfo(dcm_obj);
%---get the figures
fig_num=gca;
fig_list=findobj(fig_num,'Type','line');
[~,Sxy]=my_get_axis;
%---find cursor's position, time and peak value 

for j=1:length(cursors)
[~,sg_num(j)]=ismember(cursors(j).Target,fig_list);
sg_num(j)=length(x_sg)+1-sg_num(j);
x(j)=x_sg(sg_num(j));
t(j)=cursors(j).Position(1);
Sxy_p(j)=cursors(j).Position(2);
Sxy_i(j)=mean(Sxy(1:100,sg_num(j)));
end

%------sort the data by x_sg
[x,s]=sort(x);
sg_num=sg_num(s);
t=t(s);
Sxy_p=Sxy_p(s);
Sxy_i=Sxy_i(s);

f=struct('t',t,'sg_num',sg_num,'x',x,'Sxy_p',Sxy_p,'Sxy_i',Sxy_i);

