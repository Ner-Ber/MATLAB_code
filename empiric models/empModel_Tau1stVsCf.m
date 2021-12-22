function Tau1stVsCf = empModel_Tau1stVsCf()
    
    
    %--***NOTICE!!! this is tau1st=f(Cf)***
    Tau1stVsCf = struct;
    Tau1stVsCf.form = 'pp';
    Tau1stVsCf.breaks= [0,611.865068053674,1223.73013610735];
    Tau1stVsCf.coefs = [-0.00402163871821561,7.22793173254825,-4359.63185295701,1400000;-0.00148552145273858,-0.154169011476604,-31.4435442329665,517244.255311149];
    Tau1stVsCf.pieces = 2;
    Tau1stVsCf.order = 4;
    Tau1stVsCf.dim=1;
    
end