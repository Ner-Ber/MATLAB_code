function N=anls_calcNFromSyyA(calI,calSyy,exp)
%This function gets calibration curves. applies them to each line and
%position seperetly and then integrates the total Normal Force

p=0.98;
lineNum=1:8;
Tstart=104;
Tend=206;

for j=1:length(lineNum)
    ph=phantomGetLines(exp,0,Tstart,'min',Tend,1,1,'all',lineNum(j));
    
    if(exist([exp '\Ph\a.tif']))
        a=double(imread([exp '\Ph\a.tif']));
        line0=a(lineNum(j),:);
        ph.lines=ph.lines./repmat(line0,length(ph.lines(:,1)),1);
    end
%     plot(ph.t,mean(ph.lines,2));
%     hold all
%     
    index=logical((ph.x>0).*(ph.x<150));
    Syy=csaps(calI,calSyy,p,ph.lines);
    N(:,j)=trapz(ph.x(index),Syy(:,index)')'*1E-3*1E6/9.8;
end

N=trapz(lineNum,N')'*(6E-3)/8;

