function [frontRaw front]=mixVK_calc_frontXT_from_XAN(phE,tStart,tEnd,smtX,xStart,xEnd)
%The function finds the location in each time step where A drops under the treshold.
%then normlize by the drop, A(x=xf), and find bu csap x1,x2 when A(x=x1)=A1
%and A(x=x2)=A2. that way x_c can be calculated
%Xstart,Xend in mm important that the data is step like function. Corners
%of the block should be excluded

x=phE.x;
t=phE.t;
A=phE.lines;
dx=mean(diff(x(200:600)));
xf=22;%[mm] To deffine how far from front location to take the residual value of A.

%--norlmize A and cut the data
firstLine=repmat(mean(A(1:30,:),1),length(A(:,1)),1);
Atmp=A;
indexA(A==0)=1; indexA(A~=0)=0;
A(Atmp~=0)=A(Atmp~=0)./firstLine(Atmp~=0);
A(Atmp==0)=1;
% A=A./firstLine;

[~,pixStart]=min(abs(phE.x-xStart));
[~,pixEnd]=min(abs(phE.x-xEnd));
[~,t_indexStart]=min(abs(phE.t-tStart));
[~,t_indexEnd]=min(abs(phE.t-tEnd));

x=x(pixStart:pixEnd);
t=t(t_indexStart:t_indexEnd);
A=A(t_indexStart:t_indexEnd,pixStart:pixEnd);

[A,x]=my_smooth(A',smtX,x);
A=A';

%Exceptions


% for i=1:length(t)
%    A(i,A(i,:)>1.02)=1;
%    A(i,A(i,:)<0)=0;
% end

A0=-0.07;%boundery for interp
A1=-0.12;%A value for x1
A2=-0.6;%A value for x2
A3=-0.7;%boundery for interp

%-----create locations for test
%x_test=linspace(60,130,10);
x_test=(40:10:130);
k=0;
l=0;
m=1;
for j=1:length(t)
    
    %     if(j==1)
    %         figure;
    %     end
    %----find x where A crosses threshold
    
    tresh_loc_tmp=find(A(j,:)>0.945 & A(j,:)<0.97);
    
    if ~isempty(tresh_loc_tmp)&&(tresh_loc_tmp(1)>20)
        
        k=k+1;
        %-------Front location by raw data
        tresh_loc(k)=tresh_loc_tmp(1);
        t_index(k)=j;
        
        stage1.t(k)=t(j);
        stage1.A(k)=A(j,tresh_loc(k));
        stage1.x(k)=x(tresh_loc(k));
        
        %find A_f - the drop of A
        %-----Front location after interpolation of A(:,j)
        [~,index_f1]=min(abs(x-stage1.x(k)--1.25*xf));
        
        if(index_f1>10) % only locations where contact drop is resolved   
            
            l=l+1;
            [~,index_f2]=min(abs(x-stage1.x(k)--0.8*xf));
            stage2.Af(l)=mean(A(j,index_f1:index_f2));
            stage2.t(l)=stage1.t(k);
            stage2.xRaw(l)=stage1.x(k);
            
            An=(A(j,:)-1)/(1-stage2.Af(l));
            
            [~,index0]=min(abs(An-A0));
            [~,index3]=min(abs(An-A3));
            xx=linspace(x(index3),x(index0),1000);
            yy=csaps(x(index3:index0),An(index3:index0),0.999,xx);

             
            [stage2.A1(l),index_tmp]=min(abs(yy-A1));
            stage2.x1(l)=xx(index_tmp);
            
            [~,index2]=min(abs(An-A2));
            stage2.x2Raw(l)=x(index2);
            
            [stage2.A2(l),index_tmp]=min(abs(yy-A2));
            stage2.x2(l)=xx(index_tmp);
            
            if(m<=length(x_test))
                if(stage2.xRaw(l)>x_test(m))
                    m=m+1;
                end
            end
            
        end
    end
    
end

frontRaw=stage1;
frontRaw.note=['smtX=' num2str(smtX) 'mm' ',tresh:0.95<A<0.97 '];
front=stage2;
front.note=['A1=' num2str(A1) ',A2=' num2str(A2) ',xf=' num2str(xf) ',smtX=' num2str(smtX) 'mm'];

figure;
imagesc(phE.x,phE.t,subtruct_norm(phE.lines,1));
caxis([0.7 1.1]);
ylim([t(1) t(end)]);
set(gca,'YDir','normal');
hold on;

plot(frontRaw.x,frontRaw.t,'black','LineWidth',2.5);
plot(front.xRaw,front.t,'red','LineWidth',2.5);






