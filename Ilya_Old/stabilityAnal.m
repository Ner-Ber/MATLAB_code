clear
a1=1-logspace(-3,-1);
a2=logspace(0,-4);
a=[a1 a2(2:end)];


f=@(x) a(1)*x-log(1+x);
t(1)=fzero(f,1E-2);

for j=2:length(a)
    f=@(x) a(j)*x-log(1+x);
    t(j)=fzero(f,[t(j-1) 1E6]);
end

f=1;
figure(f);
loglog(a,t,'.-');
xlabel( '\tau / \beta K / N  V');
ylabel('t / \tau');

f=f+1;
figure(f);
loglog(a,a.*t,'.-');
xlabel( '\tau / \beta K / N  V');
ylabel('\Delta \mu / \beta');

f=f+1;
figure(f);
loglog(a./dmu,dmu,'.-');
xlabel( '\tau / d_c/ \beta  V');
ylabel('(K d_c) / (N \beta)');

f=f+1;
figure(f);
loglog(dmu.^-1,a./dmu,'.-')
xlabel('(N \beta) / (K d_c)');
ylabel( '\tau / d_c/ \beta  V');

f=f+1;
figure(f);
loglog(dmu.^-1,a,'.-')
xlabel('(N \beta)/ (K d_c)');
ylabel( '\tau K / ( \beta N )   V');