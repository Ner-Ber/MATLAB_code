function Tau2ndKinkVsCf = empModel_Tau2ndKinkVsCf()
    
    
    %--***NOTICE!!! this is tau2nd=f(Cf)***
%     Tau2ndKinkVsCf = struct;
%     Tau2ndKinkVsCf.form = 'pp';
%     Tau2ndKinkVsCf.breaks= [0.0166000000000000,0.520000000000000,0.673000000000000];
%     Tau2ndKinkVsCf.coefs = [-12941081.0270547,13989301.3261608,-6253712.65781314,2000000.00000000;5986513.20677241,-5554319.24089720,-2007542.67609147,746071.928447256];
%     Tau2ndKinkVsCf.pieces = 2;
%     Tau2ndKinkVsCf.order = 4;
%     Tau2ndKinkVsCf.dim=1;

    Tau2ndKinkVsCf = struct;
    Tau2ndKinkVsCf.form = 'pp';
    Tau2ndKinkVsCf.breaks= [0,1020,1225];
    Tau2ndKinkVsCf.coefs = [-0.000333978079400719,1.17787507696142,-1980.24266108292,2000000.00000000;-0.0485975647370100,0.155902153995255,-619.789885507147,851193.506081403];
    Tau2ndKinkVsCf.pieces = 2;
    Tau2ndKinkVsCf.order = 4;
    Tau2ndKinkVsCf.dim=1;

    
end