function [acqS phS A]=anlsAcalibration(exp,Tstart,Tinterval,Tend,lineNum,Asmt,acqSsmt)

acqS=acq132_stream_load_mat([exp '\stream'],Tstart,'min',Tend,acqSsmt,'Sxx','Syy','Sxy','Uxx','Uyy','Uxy','N','F','x_sg');

phS=phantomGetLines(exp,0,Tstart,Tinterval,Tend,1,1,'all',lineNum);
A=get_A_at_x(phS,acqS.x_sg,Asmt);
acqS=synchronizeSlowPhAcq(acqS,A.t);
