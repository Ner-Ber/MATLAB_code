%% investigate sigmoid cohesive zone function

sigmoidVar = @(x,n1,n2) (1 - x.^2).*(0.5*tanh(n1*x + n2) + 0.5) + (x + 1).^2.*(0.5*tanh(-n1*x - n2) + 0.5);
Xmodel = -1:1e-3:0;
figure;
%% change n1
% N2 = 1;
% N1 = 0:1:15;
% nn = length(N1);
% Colors = MyVaryColor(nn);
% for i=1:nn
%     thisSig = @(x) sigmoidVar(x,N1(i),N2);
%     CohesiveModelStruct_i = CrackSolutionGeneralCohesive_InsertVariables_Neri(...
%         v,PhediStruct_7_12.Param.CohsvModl_h,PhediStruct_7_12.solAtSG.Gamma,'L',10e-3,thisSig);
%     subplot(2,2,1); hold on;
%     plot(Xmodel,thisSig(Xmodel),'Color',Colors(i,:),'LineWidth',1.5,'DisplayName',['n1=',num2str(N1(i)),'  n2=',num2str(N2)]);
%     subplot(2,2,2); hold on;
%     plot(CohesiveModelStruct_i.x(2:end)+0.003,diff(CohesiveModelStruct_i.Ux)./diff(CohesiveModelStruct_i.t),...
%         'Color',Colors(i,:),'LineWidth',1.5,'DisplayName',['n1=',num2str(N1(i)),'  n2=',num2str(N2)]);
%     pause(0.1);
% end
% 
%% change n2
% N2 = -1:0.8:10;
% N1 = 8;
% nn = length(N2);
% Colors = MyVaryColor(nn,parula);
% for i=1:nn
%     thisSig = @(x) sigmoidVar(x,N1,N2(i));
%     CohesiveModelStruct_i = CrackSolutionGeneralCohesive_InsertVariables_Neri(...
%         v,PhediStruct_7_12.Param.CohsvModl_h,PhediStruct_7_12.solAtSG.Gamma,'L',10e-3,thisSig);
%     subplot(2,2,3); hold on;
%     plot(Xmodel,thisSig(Xmodel),'Color',Colors(i,:),'LineWidth',1.5,'DisplayName',['n1=',num2str(N1),'  n2=',num2str(N2(i))]);
%     subplot(2,2,4); hold on;
%     plot(CohesiveModelStruct_i.x(2:end)+0.003,diff(CohesiveModelStruct_i.Ux)./diff(CohesiveModelStruct_i.t),...
%         'Color',Colors(i,:),'LineWidth',1.5,'DisplayName',['n1=',num2str(N1),'  n2=',num2str(N2(i))]);
%     pause(0.1);
% end

%% 
n1Vec = 5:2:35;
L_vec = 3:2:33;
for n_i = 1:length(n1Vec)
    for L_i = 1:length(L_vec)
        CohesiveModelStructIntrfc = CrackSolutionGeneralCohesive_InsertVariables_Neri(PhediData.Cf/DataStruct.Param.Cr,1e-7,DataStruct.solAtSG.Gamma,'L',6e-3,cohsvMdl);
        
        
    end
end