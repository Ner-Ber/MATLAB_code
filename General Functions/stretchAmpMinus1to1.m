function strechedMeasureMatrix = stretchAmpMinus1to1(measureMatrix)
%sterchedMeasureMatrix = stretchAmpMinus1to1(measureMatrix)
%
%measureMatrix is a mtrix where each column is a measurement

measureMatrixZeroed = bsxfun(@minus,measureMatrix, min(measureMatrix));
strechedMeasureMatrix = bsxfun(@rdivide,measureMatrixZeroed, max(measureMatrixZeroed))*2-1;
end