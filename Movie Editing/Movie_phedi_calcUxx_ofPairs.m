function [UxxFromPhedi,phediPairs,initial_locations,L_vec,u_mat] = Movie_phedi_calcUxx_ofPairs(PhediLocation, varargin)
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


%% set defaults
[initialPhediLocations,distThresh] =...
    setDefaults4function(varargin,PhediLocation(1,:),150e-6);


%% select phedis by slope incline:
% relevantPhediLocations = initialPhediLocations(slope_incline_vec==side);
relevantPhediLocations = initialPhediLocations;

allPairs = [(2:length(initialPhediLocations))',(1:length(initialPhediLocations)-1)'];  % all pairs of phedis
%--- select only pairs of the same asperity:
pairsInitialVals = relevantPhediLocations(allPairs);
pairsDist = diff(pairsInitialVals,1,2);
pairsLogical = abs(pairsDist)<distThresh;
selectedPairs = allPairs(pairsLogical,:);
%--- measure strain in all pairs of phedis
m = size(selectedPairs,1);
L_vec = [];
u_mat = [];
epsilon_mat = [];
initial_locations = [];
phediPairs = [];
for pairIdx = 1:m
    ValPositive = relevantPhediLocations(selectedPairs(pairIdx,1));
    ValNegative = relevantPhediLocations(selectedPairs(pairIdx,2));
    L = ValPositive-ValNegative;
%     %--- set relevant distances to consider
%     if ~(minDistance4calc<L && L<maxDistance4calc)||...
%             ~((slope_incline_vec(selectedPairs(pairIdx,1))==side)&&(slope_incline_vec(selectedPairs(pairIdx,2))==side))
%         continue
%     end
    
%     L_prime = PhediLocation(:,j)-PhediLocation(:,i);
    L_prime = PhediLocation(:,selectedPairs(pairIdx,1))-PhediLocation(:,selectedPairs(pairIdx,2));
    u_x =  L_prime-L;
    epsilon = u_x/L;
    
    L_vec = cat(2,L_vec,L);
    u_mat = cat(2,u_mat,u_x);
    epsilon_mat= cat(2,epsilon_mat,epsilon);
    initial_locations = cat(2,initial_locations,[ValNegative; ValPositive]);
    phediPairs = cat(2,phediPairs,selectedPairs(pairIdx,:)');
end

UxxFromPhedi = epsilon_mat;

end