function Ux = Crack_particleDisplacement(x,vx,Cf)
%Ux = Crack_particleDisplacement(x,vx)
%
% Crack_particleDisplacement will return a particle displacement, Ux,
% given the particle velocity, location, and front velocity
% this basically means integrating over x. 

Ux = -(1/Cf)*cumtrapz(x,vx);

end