function f=get_value_at_x_sg(f_in,x_in,x_sg,smtX)
%The function gets function f_in and coordinate, smoothes the data and
%extract value at x_sg position

if nargin==2
    smt=1;
end


for sg_num=1:length(x_sg)
    [~,index(sg_num)]=min(abs(x_in-x_sg(sg_num)));

end
f_in=smooth(f_in,smtX);
f=f_in(index);
   


