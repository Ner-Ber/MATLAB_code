function [X,GammaVec] = gamma_calcGammaAlongInterface(DataStruct)
%gamma_calcGammaAlongInterface will calculate the fracture energy from all
%the strain gages along the interface

sgDataStruct = DataStruct.SgData;
BigPicRotStruct = DataStruct.BigPicRotStruct;
N = length(sgDataStruct.x_sg);
relevantSG = 1:N;
GammaVec = zeros(1,N);

for i = 1:N
%     solAtSG = gamma_findGammaFit(sgDataStruct,BigPicRotStruct, relevantSG(i));
    solAtSG = gamma_findGammaFit_ByAmpUxx(sgDataStruct,BigPicRotStruct, relevantSG(i));
    GammaVec(i) = solAtSG.Gamma;
end
X = sgDataStruct.x_sg;
