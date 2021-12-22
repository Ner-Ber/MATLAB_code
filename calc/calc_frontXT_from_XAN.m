function [frontRaw front]=calc_frontXT_from_XAN(phE,tStart,tEnd,smtX,xStart,xEnd)
%The function finds the location in each time step where A drops under the treshold.
%then normlize by the drop, A(x=xf), and find bu csap x1,x2 when A(x=x1)=A1
%and A(x=x2)=A2. that way x_c can be calculated
%Xstart,Xend in mm important that the data is step like function. Corners
%of the block should be excluded
%For negative direction use xStart and xEnd from nucleation and end
%respectively.

x=phE.x;
t=phE.t;
A=phE.lines;
dx=mean(diff(x(200:600)));

%--norlmize A and cut the data
%firstLine=repmat(mean(A(1:30,:),1),length(A(:,1)),1);
[~,t_indexStart]=min(abs(phE.t-tStart));
firstLine=repmat(mean(A(t_indexStart-10:t_indexStart,:),1),length(A(:,1)),1);
A=A./firstLine;

[~,pixStart]=min(abs(phE.x-xStart));
[~,pixEnd]=min(abs(phE.x-xEnd));
[~,t_indexStart]=min(abs(phE.t-tStart));
[~,t_indexEnd]=min(abs(phE.t-tEnd));

if(pixStart>pixEnd)
    Dindex=-1; %for ruptures propagating in the negative direction
else
    Dindex=1;
end

xf= Dindex*10;%[mm] To deffine how far from front location to take the residual value of A.

t=t(t_indexStart:t_indexEnd);
x=x(pixStart:Dindex:pixEnd);
A=A(t_indexStart:t_indexEnd,pixStart:Dindex:pixEnd);

[A,x]=my_smooth(A',smtX,x);
A=A';

A0=-0.07;%boundery for interp
A1=-0.12;%A value for x1
A2=-0.63;%A value for x2
A3=-0.7;%boundery for interp

%-----create locations for test
%x_test=linspace(60,130,10);
x_test=(135:1:165);
k=0;
l=0;
m=1;
for j=1:length(t)
    
    %     if(j==1)
    %         figure;
    %     end
    %----find x where A crosses threshold
    
    tresh_loc_tmp=find(A(j,:)>0.95 & A(j,:)<0.97);
    %if ~isempty(tresh_loc_tmp)&&(tresh_loc_tmp(1)>20)
    if ~isempty(tresh_loc_tmp)&&(tresh_loc_tmp(1)>5)
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
            
            %------fit to exp
            %              [~,index1]=min(abs(An-A1));
            %             xfit=-(x(index_f1:index2)-stage1.x1(l));
            %             Anfit=An(index_f1:index2)+1;
            %
            %             [xData, yData] = prepareCurveData( xfit, Anfit );
            %
            %             % Set up fittype and options.
            %             ft = fittype( 'a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y' );
            %             opts = fitoptions( ft );
            %             opts.Display = 'Off';
            %             opts.Lower = [0 0.1 0];
            %             opts.StartPoint = [0.5 3 0.8];
            %             opts.Upper = [10 10 10];
            %
            % Fit model to data.
            %             [fitresult,~] = fit( xData, yData, ft, opts );
            %             stage1.xc_fit(l)=fitresult.b;
            %---check fit
            %             figure(13);
            %             plot(xfit,Anfit,'.-');
            %             hold all
            %             plot(fitresult);
            %             xlim([-3 20]);
            %             ylim([-0.1 1.1]);
            %             hold off;
            %------plot for check
            
            if(m<=length(x_test))
                if(stage2.xRaw(l)>x_test(m))
                    m=m+1;
                    %
                    %                                          figure(14);
                    %                                          plot(x-stage2.x2(l),(A(j,:)-1)/(1-stage2.Af(l)),'.-');
                    %                                          hold all;
                    %figure(14);
                    %                                         plot(x-stage2.xRaw(l),An,'.-');
                    %                                         hold all;
                    %figure(15);
                    %plot((x-stage2.x1(l))/(stage2.x1(l)-stage2.x2(l)),An,'.-');
                    %hold all;
                    %                                         plot(xx-stage1.x1(k),yy);
                    %                                         plot(stage1.x2(k)-stage1.x1(k),A2,'o');
                    %                                         plot([stage1.xraw(k) stage1.x1(k) stage1.x2(k)]-stage1.x1(k),[A1 A1 A2],'o');
                    
                end
            end
            
        end
    end
    
end

% figure(13);
% legend(num2str(x_test'));legend off;
% figure(14);
% legend(num2str(x_test'));legend off;

frontRaw=stage1;
frontRaw.note=['smtX=' num2str(smtX) 'mm' ',tresh:0.95<A<0.97 '];
front=stage2;
front.xf=xf;
front.note=['A1=' num2str(A1) ',A2=' num2str(A2) ',xf=' num2str(xf) ',smtX=' num2str(smtX) 'mm'];


% figure;
% imagesc(phE.x,phE.t,subtruct_norm(phE.lines,1));
% caxis([0.7 1.1]);
% ylim([t(1) t(end)]);
% set(gca,'YDir','normal');
% hold on;

% plot(frontRaw.x,frontRaw.t,'black','LineWidth',2.5);
% plot(front.xRaw,front.t,'red','LineWidth',2.5);



% front.tRaw=stage2.t;
% front.xRaw=stage2.x;
% front.diff_A_PeakRaw=stage2.diff_A_Peak;
% %---Exclude locations with large derivitive
%
% vRaw=diff(front.xRaw)./diff(front.tRaw);
% include_index=vRaw>-1000;
% include_index=include_index.*(vRaw<3500);
% include_index=logical(include_index);
% include_index=[include_index(1) include_index];
% front.t=front.tRaw(include_index);
% front.x=front.xRaw(include_index);
% front.diff_A_Peak=front.diff_A_PeakRaw(include_index);





