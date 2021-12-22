%isPrecursor = phantom_isPrecursor(BigPicRotStruct)
%
% input maybe either BigPicRotStruct or just the 2D DataMatNrom matrix. 

function isPrecursor = phantom_isPrecursor(BigPicRotStruct)
    
if isstruct(BigPicRotStruct)
    DataMatNorm = BigPicRotStruct.DataMatNorm;
else
    DataMatNorm = BigPicRotStruct;
end


%-- take central 5% of frame:
L = size(DataMatNorm,2);
croppedMat = DataMatNorm(:,round([L/2-L*0.025, L/2+L*0.025]));
centralVal = mean(croppedMat(:));


Thresh = 0.94;
isPrecursor = centralVal>=Thresh;

end