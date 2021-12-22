function RubberBand(plot_num)
%Enter the plot number as argument (i.e 1 for the first plot at the browser)

% to run this function type:
% [x_both_relevant,y_both_relevant,weights_both_relevant]=rubber_band(x_both,y_both,q,weights_both)

if nargin<1
    plot_num=1;
end
ind=[];
leng_end=0;

h=gcf;
% b=get(h,'Children');
% c=get(b,'Children');
c=get(gca,'Children');
c(1:end)=c(end:-1:1);
x_both=get(c(plot_num),'XData');
y_both=get(c(plot_num),'YData');        

% hold off
% plot(x_both,y_both,'.');
% ylim([(min(y_both)-5) (max(y_both)+5)]);
flag=0;
iter=0;
while flag ==0
    k = waitforbuttonpress;
  
    if k==0 %mouse button click
          SelectionType=get(h,'SelectionType');
        if ~strcmp(SelectionType, 'extend')
            continue;
        end
        iter=iter+1;
        point1 = get(gca,'CurrentPoint');   % button down detected
        finalRect = rbbox;
        point2 = get(gca,'CurrentPoint');    % button up detected
        point1 = point1(1,1:2);              % extract x and y
        point2 = point2(1,1:2);
%         p1 = min(point1,point2);             % calculate locations
%         offset = abs(point1-point2);         % and dimensions
%         x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
%         y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
        temp=find((x_both>min(point1(1),point2(1))) & (x_both<max(point1(1),point2(1))) & (y_both>min(point2(2),point1(2))) & (y_both<max(point2(2),point1(2))) );
        %temp=find(x_both>point1(1) & x_both<point2(1) & y_both>point2(2) & y_both<point1(2) );
        leng_start=leng_end+1;
        leng_end=leng_end+length(temp);
        if leng_start<=leng_end
            ind(leng_start:leng_end)=temp;
        end
        hold all
        %if findobj(c(plot_num),'Color','r')
            Color='b';
        %else
        %    Color='r';
        %end
            
        plot(x_both(temp),y_both(temp),[Color '.'],'DisplayName',sprintf('Part %d',iter))

        clear x y point1 point2 p1 offset finalRect temp
        

    else
        flag=1;
    end

end

clear leng_end leng_start k flag
% ind=unique(ind);
% c = setdiff(index_vector,ind);
%I want to keep those
% x_both_relevant=x_both(c);
% y_both_relevant=y_both(c);
x_both_relevant=x_both(ind);
y_both_relevant=y_both(ind);

plot(x_both_relevant,y_both_relevant,'k.','DisplayName','Rubbered');

%grid minor








