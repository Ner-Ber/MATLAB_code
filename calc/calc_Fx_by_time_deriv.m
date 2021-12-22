function spl=calc_Fx_by_time_deriv(struct,frontV,SxxSmt)
%frontV is row oriented front velocity vector 
%SxxSmt should be avector
% for slow front good result is about smt=151 for Cf=157; 
if (nargin<3)
SxxSmt=floor(151./frontV*157);
end

%------------Read exp_details.txt
expDetails=expDetailsRead(struct.exp);
height=expDetails.sg_height;%[mm]
%---------------
Fx0=calc_Fx_stress_splined(struct);
spl.Fx0=mean(Fx0.Fx(1:40,2:2:end-1),1);
clear Fx0;


Sxy=my_smooth(subtruct_norm(struct.Sxy),5);
%Sxx=subtruct_norm(struct.Sxx);
%Sxx=my_smooth(Sxx,SxxSmt);
Uxx=subtruct_norm(struct.Uxx);
Sxx=my_smooth(-3*Uxx,SxxSmt);
t=struct.t;
clear struct;

fudgeFactor=1;
spl.t=t(2:end);
dSxxdx=diff(Sxx,1,1)./repmat(diff(t),1,length(Sxx(1,:)));
dSxxdx=--dSxxdx./repmat(frontV,length(t)-1,1)*1000;%'-' is needed because positive Sxx is compresion. [frontV]=[m/s] [t]=[msc] -> 1000 factor needed
%another '-' is needed because dx is - dt
spl.dSxxdx=dSxxdx.*repmat(height,length(t)-1,1)*fudgeFactor/1000; %1000 convert from [h]=[mm] -> to m
spl.Sxy=Sxy(2:end,:);
spl.Sxx=Sxx(2:end,:);
spl.Fx=spl.Sxy+spl.dSxxdx;
spl.Fx0=repmat(spl.Fx0,length(spl.t),1);
