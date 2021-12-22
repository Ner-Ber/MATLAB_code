function NuVsFreq = empModel_NuVsFreq()
    
    
    %--***NOTICE!!! this is Nu=f(log10(freq))***
    %-- based on data from Read&Duncan 1981 Fig5 --
    %-- here freq is given in Hz --%
    
    NuVsFreq = struct;
    NuVsFreq.form = 'pp';
    NuVsFreq.breaks= [-2.11757269114235,1.66245259244787,5.44247787603809];
    NuVsFreq.coefs = [0.000674019957070780,-0.00486556350909626,-2.67371449937667e-25,0.367979798471214;-0.000305823199385989,0.00277787392901958,-0.00789151939697770,0.334862401632551];
    NuVsFreq.pieces = 2;
    NuVsFreq.order = 4;
    NuVsFreq.dim=1;
    
end