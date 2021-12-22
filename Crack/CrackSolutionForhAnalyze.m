function anl=CrackSolutionForhAnalyze(v,h)
%h[m]

%-------Materila parameters
sol=CrackSolutionGeneralCohesive(0.01,h);
anl.Cr=sol.Cr;
anl.Cs=sol.Cs;
anl.Cd=sol.Cd;
anl.Gamma=sol.Gamma;
anl.tau_p=sol.tau_p;
anl.y=sol.y;

%------Define V
%-----Sub-Rayleigh - v[Cr]
% v=[linspace(0.001,0.8,20) linspace(0.81,0.96,15) linspace(0.962,0.997,10)];
% anl.v=v*anl.Cr;

%---Supershear  - v[m/s]
v=linspace(anl.Cs*1.4,anl.Cd*0.995,60);
anl.v=v;


%solTmp=Crack_EqMotion(v);
%anl.l=solTmp.l;
Uyymax=zeros(length(v),1);
Uxxmax=Uyymax;
Xc=Uyymax;


parfor j=1:length(v)
    %waitbar(j / length(v))
    %---------Analyze Slip weakening model
    %sol=CrackSolutionGeneralCohesive(v(j),h);

    %-----------Analyze LEFM
    %sol=CrackSolutionForh(v(j),-1/2,h);
    %sol=CrackSolutionSmoothedSg(v(j));
    %sol=CrackSolutionSmoothedByMaskSg(v(j),h);
    %-----Analyze Broberg
    %     sol=CrackSolution_SelfSimil arFreund(v(j),51e-3,-3.5e-3);
    %     anl.Uxy0(j)=sol.Uxy0;
    %     anl.l(j)=sol.l;
    
    %-------------SuperShear
     sol=Crack_SupershearCohesiveForhNew(v(j),h);
    
    
    %--------------Uxy
    %sol.Uxy=smooth(sol.Uxy,sol.smtT);
    % [indexMax indexMin]=localMaxMin(sol.Uxy);
    %     if ~isempty(indexMax)
    %         anl.Uxy.x_max(j,:)=sol.x(indexMax(end));
    %         anl.Uxy.max(j,:)=sol.Uxy(indexMax(end));
    %     end
    %     if ~isempty(indexMin)
    %         anl.Uxy.x_min(j,:)=sol.x(indexMin);
    %         anl.Uxy.min(j,:)=sol.Uxy(indexMin);
    %     end
    %anl.Uxy.max(j,:)=max(abs(sol.Uxy));
    
    %--------------Uxy
    %sol.Uxy=smooth(sol.Uxy,sol.smtT);
    
    
%     [indexMax indexMin]=localMaxMin(sol.Uxy);
%     if ~isempty(indexMax)
%         anl.Uxy.x_max(j,:)=sol.x(indexMax);
%         anl.Uxy.max(j,:)=sol.Uxy(indexMax);
%     end
%     if ~isempty(indexMin)
%         anl.Uxy.x_min(j,:)=sol.x(indexMin);
%         anl.Uxy.min(j,:)=sol.Uxy(indexMin);
%     end
% %     
%             if (j==1)
%                 figure;
%                 plot(sol.x,sol.Uxy,'.-');
%                 hold all;
%                 for k=1:length(anl.Uxy.x_max)
%                     plot(anl.Uxy.x_max(j,k),anl.Uxy.max(j,k),'o');
%                 end
%                 for k=1:length(anl.Uxy.x_min)
%                     plot(anl.Uxy.x_min(j,k),anl.Uxy.min(j,k),'o');
%                 end
%             end
%     
    
    %--------------Uyy
   %sol.Uyy=smooth(sol.Uyy,sol.smtT);
   
   
%     [indexMax indexMin]=localMaxMin(sol.Uyy);
%     if ~isempty(indexMax)
%         anl.Uyy.x_max(j,:)=sol.x(indexMax);
%         anl.Uyy.max(j,:)=sol.Uyy(indexMax);
%     end
%     if ~isempty(indexMin)
%         anl.Uyy.x_min(j,:)=sol.x(indexMin);
%         anl.Uyy.min(j,:)=sol.Uyy(indexMin);
%     end
%     
    
Uyymax(j,:)=max(abs(sol.Uyy));
 %Uyymax(j,:)=sol.Uyy;   
     
    %      %--------------Uxx
    %sol.Uxx=smooth(sol.Uxx,sol.smtT);
    
%     [indexMax indexMin]=localMaxMin(sol.Uxx);
%     
%     if ~isempty(indexMax)
%         anl.Uxx.x_max(j,:)=sol.x(indexMax);
%         anl.Uxx.max(j,:)=sol.Uxx(indexMax);
%     end
%     if ~isempty(indexMin)
%         anl.Uxx.x_min(j,:)=sol.x(indexMin);
%         anl.Uxx.min(j,:)=sol.Uxx(indexMin);
%     end
    Uxxmax(j,:)=max(abs(sol.Uxx));
    %Uxxmax(j,:)=sol.Uxx;
    %------Xc
    
    Xc(j)=sol.Xc;
    
end
anl.Uyy.max=Uyymax;
anl.Uxx.max=Uxxmax;
anl.Xc=Xc;




%anl.Xc0=sol.Xc0;

%[num2str(7.5) ',' num2str(anlS.Cd) ',' num2str(anlS.Cs) ',' num2str(anlS.Gamma) ',' num2str(anlS.tau_s/1e6) ',' num2str(anlS.Xc0*1e3) ]

function [indexMax indexMin]=localMaxMin(f)

%--find Maximum
indexMax=localMaximum(f,10);
tmpindex=indexMax>10;
indexMax=indexMax(tmpindex);
tmpindex=indexMax<length(f)-10;
indexMax=indexMax(tmpindex);

%--find Min
indexMin=localMaximum(-f,10);
tmpindex=indexMin>10;
indexMin=indexMin(tmpindex);
tmpindex=indexMin<length(f)-10;
indexMin=indexMin(tmpindex);
