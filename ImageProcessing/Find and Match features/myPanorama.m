
warning('OFF','images:initSize:adjustingMag');


% create panorama from my parent's house!
tic;
numFrames = 4;
inpPathFormat = 'Home%d.jpg';
outPath = 'data/out/mine/Home.jpg';
renderAtFrame = ceil(numFrames/2);
generatePanorama(inpPathFormat,outPath,numFrames,renderAtFrame,true);
toc;
pause(2);
close all;

% create panorama from Mitzpe Shaharoot!
tic;
numFrames = 3;
inpPathFormat = 'Negev%d.png';
outPath = 'data/out/mine/Negev.jpg';
renderAtFrame = ceil(numFrames/2);
generatePanorama(inpPathFormat,outPath,numFrames,renderAtFrame,true);
toc;
pause(2);
close all;



warning('ON','images:initSize:adjustingMag');
