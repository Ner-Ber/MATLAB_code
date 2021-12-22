function numberOfPix=calcEffectiveArea(exp_dir,Tstart,Tend,smtX,smtT)
%smtT=5


for j=1:8
    phS(j)=phantomGetLines(exp_dir,0,Tstart,'min',Tend,1,smtX,'all',j);
    phS(j).lines=my_smooth(phS(j).lines,smtT);
    tmp(j).bin=phS(j).lines*0;
    tmp(j).bin(phS(j).lines>50)=1;
end

numberOfPix=tmp(1).bin;
for j=2:8
    numberOfPix=numberOfPix+tmp(j).bin;
end

