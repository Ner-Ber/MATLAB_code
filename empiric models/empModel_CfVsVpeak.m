function CfVsVpeak = empModel_CfVsVpeak()
    
    
    %--***NOTICE!!! this is Cf=f(Vpeak)***
    CfVsVpeak = struct;
    CfVsVpeak.form = 'pp';
    CfVsVpeak.breaks= [0.0224954422506689,0.502802386390289,0.983109330529909,1.46341627466953];
    CfVsVpeak.coefs = [-800.495594714734,-1349.69775285275,3052.82844700490,17.9562182468574;1737.18671783982,-2503.14853153674,1202.27962191009,1084.18437023562;-1.39431978095863e-10,-1.01782155462987e-12,1.03190566880141e-10,1276.67212063598];
    CfVsVpeak.pieces = 3;
    CfVsVpeak.order = 4;
    CfVsVpeak.dim=1;
    
end