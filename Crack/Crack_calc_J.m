function [theta,J]=Crack_calc_J(sol)

jx=sol.Sxy.*sol.vy+sol.Sxx.*sol.vx;
jy=sol.Sxy.*sol.vx+sol.Syy.*sol.vy;

theta=sol.theta;
r=sol.r;
J=(jx.*cos(theta)+jy.*sin(theta))*r; 
J=trapz( jx.*cos(theta)+jy.*sin(theta) )*r*(theta(2)-theta(1));