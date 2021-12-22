function EvsFreq = empModel_EvsFreq()
    
    
    %--***NOTICE!!! this is E=f(log10(freq))***
    %-- based on data from Read&Duncan 1981 Fig4 --
    %-- here freq is given in Hz, E in GPa --%
    
    EvsFreq = struct;
    EvsFreq.form = 'pp';
    EvsFreq.breaks= [-2.04847131712159,2.24381574221883,6.53610280155925];
    EvsFreq.coefs = [-0.0116373449169612,0.105811259380717,0.154681013696021,3.37882420409268;-3.82547683417129e-05,-0.0440412155957442,0.419815773289151,5.07191697780237];
    EvsFreq.pieces = 2;
    EvsFreq.order = 4;
    EvsFreq.dim=1;
    
end