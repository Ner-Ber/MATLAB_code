function spl=calc_Fx_stress_splined(struct)
%The function gets struct that conatins t,x_sg,Sxy,Sxx,Syy
%The function calculates the total force on volume element of sg, and
%splined stress.
% notice that (Fx dSxxdx)'s   x dimention is smaller because of the derivitive. the appropriate x vector is x_spl(2:end-1); 
% plot(Fx.x(3:2:end-2),Fx.Fx(:,2:2:end-1)./acqS.Syy(:,2:end-1),'.-');

%------------Read exp_details.txt
expDetails=expDetailsRead(struct.exp);
height=expDetails.sg_height;%[mm]
%---------------
x_sg=struct.x_sg;
Sxy=struct.Sxy;
Sxx=struct.Sxx;
Syy=struct.Syy;
spl.t=struct.t;
clear struct;

fudgeFactor=1;

spl.x=interp1(1:length(x_sg),x_sg,1:0.5:length(x_sg)); %x coordinate of the interpolation
spl.Sxy=spline(x_sg,Sxy,spl.x);  %interpolate the stresses 
spl.Sxx=spline(x_sg,Sxx,spl.x);
spl.Syy=spline(x_sg,Syy,spl.x);
dSxxdx=-diff(Sxx,1,2)./repmat(diff(x_sg),length(Sxx(:,1)),1)*1000; %differentiate along x axis.'-' is needed because positive Sxx is compresion.
%factor 1000 is from mm to meter of x_sg
height=(height(1:end-1)+height(2:end))/2; %mean hight betwean two sg. approximation may be improved
height=height*fudgeFactor;
dSxxdx=dSxxdx.*repmat(height,length(Sxx(:,1)),1)/1000;
%first h*1000 -> elemnt hight 

%dSxxdx=-diff(Sxx,1,2)./repmat(diff(x_sg),length(Sxx(:,1)),1)*4;

x_d_spl=spl.x(2:end-1); % x coordinate matched to dSxx_dx_spl  
spl.dSxxdx=spline(x_d_spl(1:2:end),dSxxdx,x_d_spl); %interpolate the derivative
spl.Fx=spl.Sxy(:,2:end-1)+spl.dSxxdx(:,:); %the total force (interpolated)


%calculate Syy at the interface
dSxydx=diff(Sxy,1,2)./repmat(diff(x_sg),length(Sxy(:,1)),1)*1000; %differentiate along x axis.'-' is needed because positive Sxx is compresion.
dSxydx=dSxydx.*repmat(height,length(Sxy(:,1)),1)/1000;
spl.dSxydx=spline(x_d_spl(1:2:end),dSxydx,x_d_spl); %interpolate the derivative
spl.Syy0=spl.Syy(:,2:end-1)+spl.dSxydx(:,:); %the total force (interpolated)

if nargout==0
    field_names=fieldnames(spl);
    for j=1:length(field_names)
        assignin('caller',field_names{j},spl.(field_names{j}));
    end
    clear spl;
end
