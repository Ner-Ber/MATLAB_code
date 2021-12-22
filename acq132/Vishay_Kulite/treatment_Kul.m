function K=treatment_Kul(ch_Master,ch_Slave,ch_offset_Master,ch_offset_Slave,exp_details,sg_names)

axis_sort_index=exp_details.axis_sort_index;

ch_Master=calculate_gagefactor_Kul(ch_Master,ch_offset_Master,exp_details,'Master');
ch_Slave=calculate_gagefactor_Kul(ch_Slave,ch_offset_Slave,exp_details,'Slave');

sg=[ch_Master(:,exp_details.sg_Master_order) ch_Slave(:,exp_details.sg_Slave_order)];
sgtmp=sg; sgtmp(:,~sg_names)=0;

%------order the strain gages
U1=sgtmp(:,3:3:end);
U1=U1(:,axis_sort_index)*1000; %[mStrain]
U2=sgtmp(:,2:3:end);
U2=U2(:,axis_sort_index)*1000; %[mStrain]
U3=sgtmp(:,1:3:end);
U3=U3(:,axis_sort_index)*1000; %[mStrain]
    
%[K_Sxx,K_Syy,K_Sxy,K_Uxx,K_Uyy,K_Uxy,K_note]=Kul_calculate_stress_strain(sgK,exp_detailsK);

[Sxx,Syy,Sxy,Uxx,Uyy,Uxy,note]=calculate_stress_strain(U1,U2,U3,exp_details.sg_angle);
K=struct('Uxx',Uxx,'Uyy',Uyy,'Uxy',Uxy,'Sxx',Sxx,'Syy',Syy,'Sxy',Sxy,'note',note);

function ch=calculate_gagefactor_Kul(ch,ch_offset,exp_details,type)

if strcmp(type,'Master')
    ch_names=exp_details.ch_Master_names;
    Vref=exp_details.Vref_Master;
    factor=exp_details.factor_Master;
elseif strcmp(type,'Slave')
    ch_names=exp_details.ch_Slave_names;
    Vref=exp_details.Vref_Slave;
    factor=exp_details.factor_Slave;
else
    return
end

ch_K1=strcmp(ch_names,'sgK1');
ch_K2=strcmp(ch_names,'sgK2');
ch_K=or(ch_K1,ch_K2);

Resistances=ones(1,length(ch(1,:)));
Resistances(ch_K1)=(ch_offset(ch_K1)-Vref(1))./factor(ch_K1); 
Resistances(ch_K2)=(ch_offset(ch_K2)-Vref(2))./factor(ch_K2);
ch=ch./repmat(Resistances,length(ch(:,1)),1);

MasterGain=ones(length(ch(:,1)),length(ch(1,:)));
MasterGain(:,ch_K)=nlGain(ch(:,ch_K));
ch=ch./MasterGain;

function g=nlGain(x)

a=5658;
b=132.85;
g=((-b+sqrt(b^2+4*a*x))/(2*a)./x).^(-1);
