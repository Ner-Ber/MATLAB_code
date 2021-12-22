function ransaceFit = RANSAC_linearFit(coordinates,varargin)
%RANSAC_linearFit will find a linear fit using the ransac algorithm.
%
% coordinates - an Nx2 matrix, fir column is X coordinates, second column Y
% coordinates.
% OPTIONAL INPUTS [defaults]:
% numberIter - maximal number of iterations [round(N*2/3)]
% inlierTol - inlier/outlier tolerance      [0.33*std(vertical scatter)]
% slopeLimit - maximal slope (in abs) that is accepted as a model [2]

%% set defaults
N = size(coordinates,1);
[numberIter , inlierTol , slopeLimit] = setDefaults4function(varargin, round(N*2/3), 0.5*std(coordinates(:,2)),30);

%% distance function
dist = @(x0,y0,m,k)  abs(k+m*x0-y0)./sqrt(1+m^2);

%% prepare coordinate pairs in random order
allPairs = nchoosek(1:N,2);
NumParis = nchoosek(N,2);
p = randperm(NumParis,min(NumParis,numberIter));
%--- shuffle order of pairs
% shuffleVec = randperm(size(allPairs,1));
allPairsShuffled = allPairs(p,:);
R1 = allPairsShuffled(:,1);
R2 = allPairsShuffled(:,2);

% RR = randperm(N);
% R1 = RR(1:floor(N/2));
% R2 = RR((floor(N/2)+1):2*floor(N/2));

%% iterate

iteration = 1;
numberOfInlieres = [];
m_vec = [];
k_vec = [];
while iteration<=numberIter && iteration<=numel(R1)
    %     iteration = iteration+1;
    %--- choose two random coordinates
    try
        XY1 = coordinates(R1(iteration),:);
        XY2 = coordinates(R2(iteration),:);
    catch
        disp(4);
    end
    %--- create model for these two points
    Diffs = diff([XY1;XY2],1,1);
    m = Diffs(2)/Diffs(1);
    k = -m*XY1(1)+XY1(2);
    if abs(m)>slopeLimit    %% pass on sloped larger than threshold
        numberOfInlieres = cat(1,numberOfInlieres,-inf);
        m_vec = cat(1,m_vec,nan);
        k_vec = cat(1,k_vec,nan);
        iteration = iteration+1;
        continue
    end
    
    %--- exclude chosen points in distance computation
    CoorIdx = 1:N;
    CoorLogic = ~(CoorIdx==R1(iteration) | CoorIdx==R2(iteration));
    
    %--- compute distance to points
    D = dist(coordinates(:,1),coordinates(:,2),m,k);
    
    %--- count inlieres
    numberOfInlieres = cat(1,numberOfInlieres,nnz(D<=inlierTol));
    m_vec = cat(1,m_vec,m);
    k_vec = cat(1,k_vec,k);
    
    %-- propogate loop
    iteration = iteration+1;
end

%% find coordinate pair whith most inlieres
[N_inlieres, amxIdx] = max(numberOfInlieres);
% XY1 = coordinates(R1(amxIdx),:);
% XY2 = coordinates(R2(amxIdx),:);
% Diffs = diff([XY1;XY2],1,1);
% m = Diffs(2)/Diffs(1);
% k = -m*XY1(1)+XY1(2);
M = m_vec(amxIdx);
K = k_vec(amxIdx);

%% insert to outpus structure
ransaceFit = struct;
ransaceFit.slope = M;
ransaceFit.intersection = K;
ransaceFit.N_inlieres = N_inlieres;
ransaceFit.N_iterations = iteration-1;

end