function s_out=Odsv_2mu_struct(exp_dir,s,p)
%s - Enter the structure with ch_094_ex
%First it loads for example ex_dir\cal\Odscal_ch_094_ex5.mat, interpolates
%it according to s.ch_094_ex(:,5) and returns the corrcted s_out

if nargin<4
    p=0.95;
end
s_out=s;
ch_ex_num=[4 5];

for j=1:length(ch_ex_num)
cal=load([exp_dir '\cal\Odscal_ch_094_ex' int2str(ch_ex_num(j)) '.mat']);
s_out.ch_094_ex(:,ch_ex_num(j))=-csaps(cal.V,cal.d,p,s.ch_094_ex(:,ch_ex_num(j)));
s_out.ch_094_ex(:,ch_ex_num(j))=s_out.ch_094_ex(:,ch_ex_num(j))-mean(s_out.ch_094_ex(1:20,ch_ex_num(j)));
end