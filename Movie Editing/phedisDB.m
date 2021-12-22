function phedislocationInPix = phedisDB(num) 

phedis = {};
%--- date 20171015 exper 140138 event 2
phedis{1} = [1315,1266,1235,1189,1157,1111,1081,1002,953,923,876,844,797,768,721,691,637,613,535,489,457,405];
%--- date 20171024 exper 181733 event 6
phedis{2} = [610,573,548,513,487,447,427,388,365,325,304,265,243,203,181,118,82,58,22];

%% get data:
phedislocationInPix = phedis{num};

end