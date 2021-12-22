function  multipLinesROT = multipLines_createMultipLinesMat(multipLineBigPicStruct)
    names = fieldnames(multipLineBigPicStruct);
    multipLinesROT = [];
    for n = names(:)'
        multipLinesROT = cat(3,multipLinesROT,multipLineBigPicStruct.(n{1}).DataMatNorm);
    end
end