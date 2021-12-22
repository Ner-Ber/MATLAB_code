function mask=KulCrackSolutionSgMask(dx)
%dx [m] ,usually  10E-6m is fine, should be dx<1.1E-4[m]
%The function returns three masks of same size that correspond to the 3 Rosette sensors.
%Later it is used by CrackSolutionSmoothedByMaskSg.
%assumption dx=dy

if (dx>1.1E-4)
    display(['dx = ' num2str(dx) 'should be dx<1E-4']);
    exit
end
%-----vishay sg parameters EA-06-015RJ-120 
%All length scales are in [m]
sg.alpha=pi/4;
sg.l=0.82E-3;%[m] 
sg.width=0.36E-3;%approxiamtly square
sg.Up_y=0.4E-3; %y location of the upper sensor's center to the center of the rossete 
sg.Plus_y=-0.2E-3;
sg.Plus_x=0.8E-3;
sg.Minus_y=-0.2E-3;
sg.Minus_x=-0.8E-3;
sg.BoxSize=2.55E-3;% The global box (mask)size.  

%---Sensor's coordinates in pixels
sgPix.l=2*ceil(sg.l/dx/2)+1; %make the size be and odd number
sgPix.w=2*ceil(sg.width/dx/2)+1; %make the size be and odd number
sgPix.Up_y=floor(sg.Up_y/dx);
sgPix.Plus_y=floor(sg.Plus_y/dx);
sgPix.Plus_x=ceil(sg.Plus_x/dx);
sgPix.Minus_y=floor(sg.Minus_y/dx);
sgPix.Minus_x=floor(sg.Minus_x/dx);

%------Create matrice of ones for the active area of upper gage
% Box of zeros boundary is needed for easier rotation 
sgBoxSize=1.2*(sg.l^2+sg.width^2)^0.5;
sgBoxSize=2*ceil(sgBoxSize/dx/2)+1; %make the size be and odd number
sgBox=zeros(sgBoxSize);
sgBoxC_index=(sgBoxSize-1)/2+1;%center index
sgBox(sgBoxC_index-(sgPix.l-1)/2:sgBoxC_index+(sgPix.l-1)/2,sgBoxC_index-(sgPix.w-1)/2:sgBoxC_index+(sgPix.w-1)/2)=ones(sgPix.l,sgPix.w);

%----Create final mask(Need to be changed for different configuration)
maskBoxSize=floor(2.1*ceil(sg.BoxSize/dx/2))+1;%Odd size is easier
maskBox=zeros(maskBoxSize);%background of the mask is zero
maskC_index=(maskBoxSize-1)/2+1;% center index

mask.sgUp=maskBox;
mask.sgUp(maskC_index-(sgBoxSize-1)/2+sgPix.Up_y:maskC_index+(sgBoxSize-1)/2+sgPix.Up_y,maskC_index-(sgBoxSize-1)/2:maskC_index+(sgBoxSize-1)/2)=sgBox;
mask.sgUpN=sum(sum(mask.sgUp)); %Number of active pixels.Needed fot normlizing the mean 

%---Rotate sgBox and put in the correct location of the mask   
mask.sgPlus=maskBox;
mask.sgPlus(maskC_index-(sgBoxSize-1)/2+sgPix.Plus_y:maskC_index+(sgBoxSize-1)/2+sgPix.Plus_y,maskC_index-(sgBoxSize-1)/2+sgPix.Plus_x:maskC_index+(sgBoxSize-1)/2+sgPix.Plus_x)=imrotate(sgBox,45,'crop');
mask.sgPlusN=sum(sum(mask.sgPlus)); %Number of active pixels.Needed fot normlizing the mean

mask.sgMinus=maskBox;
mask.sgMinus(maskC_index-(sgBoxSize-1)/2+sgPix.Minus_y:maskC_index+(sgBoxSize-1)/2+sgPix.Minus_y,maskC_index-(sgBoxSize-1)/2+sgPix.Minus_x:maskC_index+(sgBoxSize-1)/2+sgPix.Minus_x)=imrotate(sgBox,-45,'crop');
mask.sgMinusN=sum(sum(mask.sgMinus)); %Number of active pixels.Needed fot normlizing the mean


mask.size=maskBoxSize;
figure;
imagesc(mask.sgUp+mask.sgMinus+mask.sgPlus);