%%
H = linspace(1e-9,sqrt(5e-3),200).^2;
H = [-fliplr(H),H];
for i=1:length(H)
    solSingular{i} = CrackSolutionForh_GammaChange(0.3 , -0.5 , H(i), 3.5, -0.1:1e-5:0.1);
%     sol{i} = CrackSolutionGeneralCohesive_InsertVariables_Neri(0.3,H(i),3.5,'L',3e-3,'exp',[],-0.1:1e-4:0.1);
%    disp(i); 
end
%%
SxxCell = cellfun(@(A) A.vx, solSingular, 'UniformOutput',0);
SxxMat = cell2mat(SxxCell');
SxxMatAbs = abs(SxxMat);
SxxMatAbsCut = SxxMatAbs;
SxxMatAbsCut(SxxMatAbsCut>1) = nan;
[xx,yy] = meshgrid(solSingular{1}.x,H);
figure; surf(xx,yy,SxxMatAbsCut,'EdgeColor','none');