function marked = markBetweenOnes(matrixWones)
    marked = zeros(size(matrixWones));
    for ii = 1:size(matrixWones,1)
        OnesLocation = find(matrixWones(ii,:));
        marked(ii,(OnesLocation(1)+1):(OnesLocation(end)-1)) = 1;
    end
end