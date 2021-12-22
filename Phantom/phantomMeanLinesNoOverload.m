function outLine=phantomMeanLinesNoOverload(inIm,included)


inIm=inIm.*included;
num=sum(included,1);
Sum=sum(inIm,1);
outLine=Sum./num;
outLine(isnan(outLine))=0;


end