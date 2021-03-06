%
% This script generates all of the provided example panoramas (once your 
% generatePanorama.m is ready and functional). Input images are read from 
% ../data/inp/examples and resulting panorams are then saved to 
% ../data/out/examples
%

warning('OFF','images:initSize:adjustingMag');

tic;
numFrames = 2;
inpPathFormat = 'data/inp/examples/oxford%d.jpg';
 outPath = 'data/out/examples/oxford.jpg';
renderAtFrame = ceil(numFrames/2);
generatePanorama(inpPathFormat,outPath,numFrames,renderAtFrame,true);
toc;
pause(2);
close all;

tic;
numFrames = 3;
inpPathFormat = 'data/inp/examples/backyard%d.jpg';
outPath = 'data/out/examples/backyard.jpg';
renderAtFrame = ceil(numFrames/2);
generatePanorama(inpPathFormat,outPath,numFrames,renderAtFrame,true);
toc;
pause(2);
close all;


tic;
numFrames = 4;
inpPathFormat = 'data/inp/examples/office%d.jpg';
outPath = 'data/out/examples/office.jpg';
renderAtFrame = ceil(numFrames/2);
generatePanorama(inpPathFormat,outPath,numFrames,renderAtFrame,true);
toc;
pause(2);
close all;



bonus=1;

if bonus
	numFrames = 2;
	inpPathFormat = 'data/inp/examples/parallax%d.jpg';
	outPath = 'data/out/examples/parallax.jpg';
	renderAtFrame = ceil(numFrames/2);
	generatePanorama(inpPathFormat,outPath,numFrames,renderAtFrame,true);
	pause(2);
	close all;
end

warning('ON','images:initSize:adjustingMag');
