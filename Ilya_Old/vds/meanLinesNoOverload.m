function outLine=meanLinesNoOverload(inIm,included)


inIm=inIm.*included;
num=sum(included,2);
Sum=sum(inIm,2);
outLine=Sum./num;
outLine(isnan(outLine))=0;


end