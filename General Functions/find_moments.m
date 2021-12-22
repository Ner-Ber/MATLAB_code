function [center_of_mass,var_mass,skew_mass]=find_moments(xx,yy)

center_of_mass = sum(yy.*xx)/sum(yy);
var_mass = sum(yy.*(xx-center_of_mass).^2);
skew_mass = sum(yy.*(xx-center_of_mass).^3)/var_mass^(3/2);


end
