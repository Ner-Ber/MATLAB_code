function anl=tmp_LoadData

exper={'18-07-37','18-41-46','19-07-54','19-33-50'};
l=1;

for k=1:length(exper);
    fileName=dir(['C:\Frics\Analyze\2014-01-20\' exper{k} '\*.mat']);
    for m=1:length(fileName)
        tmp=load(['C:\Frics\Analyze\2014-01-20\' exper{k} '\' fileName(m).name]);
        tmp.exper=[exper{k} '\' fileName(m).name];
        pre=anls_pre_conditions(exper{k},tmp.eventN,-2,2:7,3);
        tmp.Syy2=pre.Syy;
        anl(l)=tmp;
        clear tmp;
        l=l+1;
    end
end