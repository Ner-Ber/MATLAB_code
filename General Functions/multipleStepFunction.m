function functionHandle = multipleStepFunction(BinBorders)

BinBorders = unique(BinBorders);
functionHandle = @(x) 0;
amplitude = 1;
for i=2:length(BinBorders)
    HeaviStep = @(x) amplitude*(heaviside(x-BinBorders(i-1)).*heaviside(BinBorders(i)-x));
    functionHandle = @(x) functionHandle(x) + HeaviStep(x);
    amplitude = amplitude+1;
end