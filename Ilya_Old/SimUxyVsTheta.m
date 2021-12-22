function B=SimUxyVsTheta(xVec,yVec,Uxy)

%t=100mus
%x0=0.00;
%xBump=145.8e-3;
%xBump=134.2e-3;
%x0=0.09475;%Crack tip location
%x0=0.1182;
 %x0=0.051;
 %xBump=0.1049;%70mus
 %xBump=0.1459;%100mus
% %xBump=0.1321;%70mus
% %xBump=0.134;%Second Bump



x0=0.014;
%x0=0.004;
xBump=0.03308;%
r=xBump-x0;
%r=0.005;
k=1;
for j=1:length(yVec);
    
    y=yVec(j);
    if(y>0)
        if (r>abs(y))
            x=x0+(r^2-y^2)^0.5;
            [~, index]=min(abs(xVec-x));
            B.x(k)=xVec(index);
            B.UxyBump(k)=Uxy(j,index);
            B.theta(k)=atan(y/(x-x0));
            B.y(k)=y;
            k=k+1;
        end
    end
     
end


% %theta>pi/2
for j=length(yVec):-1:1;
    
    y=yVec(j);
    if(y>0)
        if (r>abs(y))
            x=x0-(r^2-y^2)^0.5;
            [~, index]=min(abs(xVec-x));
            B.x(k)=xVec(index);
            B.UxyBump(k)=Uxy(j,index);
            B.theta(k)=atan(y/(x-x0))+pi;
            B.y(k)=y;
            k=k+1;
        end
    end
     
end

B.x0=x0;
B.r=r;