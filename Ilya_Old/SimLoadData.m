
path='G:\Frics\Analyze\Sim\SimPlaneStress\SuperShear\2016-08-04\3d_taup1_lowerprestress_334';
cd (path)

Files=dir;
Files(1)=[];
Files(1)=[];

%----3D simulations.
j=1;
y7_5{1}.Uxx=csvread([pwd '\' Files(j).name]);
j=j+1;
y7_5{1}.Uxy=csvread([pwd '\' Files(j).name]);
j=j+1;
y7_5{1}.Uyy=csvread([pwd '\' Files(j).name]);
j=j+1;
y7_5{1}.t=csvread([pwd '\' Files(j).name]);
j=j+1;
y7_5{1}.x=csvread([pwd '\' Files(j).name]);
y7_5{1}.y=7.53e-3;
y7_5{1}.z=0;

j=j+1;
y0{1}.Sxy=csvread([pwd '\' Files(j).name]);
j=j+1;
y0{1}.SlipV=csvread([pwd '\' Files(j).name]);
j=j+1;
y0{1}.t=csvread([pwd '\' Files(j).name]);
j=j+1;
y0{1}.x=csvread([pwd '\' Files(j).name]);
y0{1}.y=0;
y0{1}.z=0;

j=j+1;
y7_5{2}.Uxx=csvread([pwd '\' Files(j).name]);
j=j+1;
y7_5{2}.Uxy=csvread([pwd '\' Files(j).name]);
j=j+1;
y7_5{2}.Uyy=csvread([pwd '\' Files(j).name]);
j=j+1;
y7_5{2}.t=csvread([pwd '\' Files(j).name]);
j=j+1;
y7_5{2}.x=csvread([pwd '\' Files(j).name]);
y7_5{2}.y=7.53e-3;
y7_5{2}.z=2.5e-3;

j=j+1;
y0{2}.Sxy=csvread([pwd '\' Files(j).name]);
j=j+1;
y0{2}.SlipV=csvread([pwd '\' Files(j).name]);
j=j+1;
y0{2}.t=csvread([pwd '\' Files(j).name]);
j=j+1;
y0{2}.x=csvread([pwd '\' Files(j).name]);
y0{2}.y=0;
y0{2}.z=2.5e-3;


save('y0','y0');
save('y7_5','y7_5');
