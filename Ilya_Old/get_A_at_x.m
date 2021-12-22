function A=get_A_at_x(vds,x_sg,smt)
%The function gets vds structure (not normilized) and gets the contact area at x_sg locations

if nargin==2
    smt=1;
end
pix=zeros(1,length(x_sg));

for sg_num=1:length(x_sg)
    [~,pix(sg_num)]=min(abs(vds.x-x_sg(sg_num)));

    A.lines(:,sg_num)=mean(vds.lines(:,pix(sg_num)-floor(smt/2):pix(sg_num)+floor(smt/2)),2);
    A.pix(sg_num)=pix(sg_num);
    A.x(sg_num)=x_sg(sg_num);

end
 %A.lines=my_smooth(vds.lines',smt)'; %my_smooth is coloumn oriented
 %A.lines=A.lines(:,pix);
   
A.t=vds.t;

