function x=phantomBuildXaxis(exp_dir,Tstart,Tinterval,Tend)

if (nargin<2)
    Tstart='start';
    Tinterval=1;
    Tend='end';
end


%--------Build x axis
PhMeta = phantomReadMeta([exp_dir '\Ph']);
if (exist([exp_dir '\exp_details.txt'],'file'))
    expDetails=expDetailsRead(exp_dir);
    [leftLim rightLim]=phantomCalcLims(exp_dir,Tstart,Tinterval,Tend);
    totalLen=expDetails.UpperBlockLength;
    x=((1:1:PhMeta.ImageWidth)-leftLim)./(rightLim-leftLim)*totalLen;
else
    error('exp_details.txt doesnt exist')
end