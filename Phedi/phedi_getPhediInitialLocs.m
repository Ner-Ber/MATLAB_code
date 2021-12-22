function phedisInitialLocations = phedi_getPhediInitialLocs(signal,varargin)
%phedisInitialLocations = phedi_getPhediInitialLocs(signal,varargin)
%
%phedi_getPhediInitialLocs will create the initial location to find the
%phedis for analysis. 
%
%INPUTS (DEFAULTS):
%1. signal (mandatory)
%2. smoothParameter (0.001) - portion of the signal length to smooth by
%3. safetyDist (30) - distance in pix from edge to eliminate phedis
%4. res (1.0036e+05) - pix/meter of the frame
%5. scratchWidth (100*1e-6) - scratch width in meters

%% defaults:
[smoothParameter, safetyDist, res, scratchWidth] = setDefaults4function(varargin,0.001,30,1.0036e+05,100*1e-6);

%% get all possible phedis
phedisInitialTrench = Movie_phedi_findTrenches(signal, smoothParameter);

%% improve
%--- eliminate phedis close to edge:
phedisInitialNoEdge = phedisInitialTrench;
phedisInitialNoEdge = phedisInitialNoEdge(abs(phedisInitialNoEdge-length(signal))>safetyDist);
phedisInitialNoEdge = phedisInitialNoEdge(phedisInitialNoEdge>safetyDist);
phedisInitialNoEdge = phedisInitialNoEdge(2:(end-1));

%--- eliminate close phedis
scratchWidth_inPix = scratchWidth*res;
PhedisDist_inPix = abs(diff(phedisInitialNoEdge));
closePairs_logical = PhedisDist_inPix<0.35*scratchWidth_inPix;
startings = find(([0;closePairs_logical]-[closePairs_logical;0])==-1);
endings = find(([0;closePairs_logical]-[closePairs_logical;0])==1);
discard = [];
for k=1:length(startings)
    discard = cat(2,discard,(startings(k)+1):endings(k)); 
end
phedisInitialLocations = phedisInitialNoEdge;
phedisInitialLocations(discard) = [];

end