function CrackSolution_plot(sol,figN)

x_offset=-sol.Xc/2*1e3;
%x_offset=0;

[~, index]=min(abs(sol.x*1e3+40));
Uxy0_Offset=sol.Uxy(index);

[~, index]=min(abs(sol.x*1e3-50));
Uyy0_Offset=sol.Uyy(index);

figure(figN(1));
hold all;
plot(sol.x*1e3+x_offset,sol.Uxy*1e3-Uxy0_Offset*1e3,'black.-');
my_legend_add( ['G=' num2str(sol.Gamma) 'Cf=' num2str(sol.v/sol.Cr) ' tauP=' num2str(sol.tau_p*1e-6) ' Xc=' num2str(sol.Xc*1e3) ' Cd=' num2str(sol.Cd) ] );

figure(figN(2));
hold all;
plot(sol.x*1e3+x_offset,sol.Uxx*1e3,'black.-');
 
figure(figN(3));
hold all;
plot(sol.x*1e3+x_offset,sol.Uyy*1e3-Uyy0_Offset*1e3,'black.-');


