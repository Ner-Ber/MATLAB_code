function VresVsCf = empModel_VresVsCf()
    
    
    %--***NOTICE!!! this is Vres=f(Cf)***
    VresVsCf = struct;
    VresVsCf.form = 'pp';
    VresVsCf.breaks= [0,629.386028151361,1258.77205630272];
    VresVsCf.coefs = [2.26858348896424e-10,-4.93106016005367e-07,0.000364700757175534,-2.29661325364740e-12;1.17210274061348e-10,-6.47615905106794e-08,1.35866800760937e-05,0.0907646313532976];
    VresVsCf.pieces = 2;
    VresVsCf.order = 4;
    VresVsCf.dim=1;
    
end