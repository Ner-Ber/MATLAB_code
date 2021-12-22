function [UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx(PhediLocation, varargin)
%[UxxFromPhedi] = Movie_phedi_calcUxx(PhediLocation, varargin)
% Movie_phedi_calcUxx will calculate the strain dL/L directly from the
% phedi measurments.
%
% ##optional arg.       (default va;.) :
% initialPhediLocations (PhediLocation(1,:))
% minDistance4calc      (-inf)
% maxDistance4calc      (inf)
% slope_incline_vec     (ones(1,size(PhediLocation,2))
% side                  (1)
%
% notice that if no 'slope_incline_vec' is inserted Movie_phedi_calcUxx
% will take all phedis into calculation.


%--- set defaults
[initialPhediLocations, minDistance4calc, maxDistance4calc, slope_incline_vec, side] =...
    setDefaults4function(varargin,...
    PhediLocation(1,:), -inf, inf, ones(1,size(PhediLocation,2)), 1);


%--- select phedis by slope incline:
% relevantPhediLocations = initialPhediLocations(slope_incline_vec==side);
relevantPhediLocations = initialPhediLocations;


N = length(relevantPhediLocations);      % number of phedis
allPairs = nchoosek(1:N,2);             % all pairs of phedis

%--- measure strain in all pairs of phedis
m = size(allPairs,1);

L_vec = [];
u_mat = [];
epsilon_mat = [];
initial_locations = [];
phediPairs = [];
for phedIdx = 1:m
    [i_val, i] = min(relevantPhediLocations(allPairs(phedIdx,:)));
    [j_val, j] = max(relevantPhediLocations(allPairs(phedIdx,:)));
    L = j_val-i_val;
    %--- set relevant distances to consider
    if ~(minDistance4calc<L && L<maxDistance4calc)||...
            ~((slope_incline_vec(allPairs(phedIdx,i))==side)&&(slope_incline_vec(allPairs(phedIdx,j))==side))
        continue
    end
    
%     L_prime = PhediLocation(:,j)-PhediLocation(:,i);
    L_prime = PhediLocation(:,allPairs(phedIdx,j))-PhediLocation(:,allPairs(phedIdx,i));
    u =  bsxfun(@minus,L_prime,L);
    epsilon = bsxfun(@rdivide,u,L);
    
    L_vec = cat(2,L_vec,L);
    u_mat = cat(2,u_mat,u);
    epsilon_mat= cat(2,epsilon_mat,epsilon);
    initial_locations = cat(2,initial_locations,[j_val; i_val]);
    phediPairs = cat(2,phediPairs,allPairs(phedIdx,:)');
end

UxxFromPhedi = epsilon_mat;

end