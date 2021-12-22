I=readmraw('example.mraw',[0 20]);

for n=1:1:10
imshow(I.Images.RawImages(:,:,n),[0 3000]);
pause(.1);
end
