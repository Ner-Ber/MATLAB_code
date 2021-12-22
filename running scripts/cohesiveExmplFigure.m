%% create cohesive example figure

Xchsv = -1:1e-3:0;
CfVec = [351 1180];
Gammas = [4 3];
for v = 1:length(CfVec)
    %% create solutions
    sol(1,v) = CrackSolutionGeneralCohesive_InsertVariables_Neri(CfVec(v)/1255,1e-8,Gammas(v),'L',3e-3,'exp');
    sol(2,v) = CrackSolutionGeneralCohesive_InsertVariables_Neri(CfVec(v)/1255,1e-8,Gammas(v),'L',3e-3,'linear');
    sol(3,v) = CrackSolutionGeneralCohesive_InsertVariables_Neri(CfVec(v)/1255,1e-8,Gammas(v),'L',3e-3,'Ohnaka89');
    sol(4,v) = CrackSolutionGeneralCohesive_InsertVariables_Neri(CfVec(v)/1255,1e-8,Gammas(v),'L',3e-3,'Const');
    HandSadParab = @(x) 1-x.^2;
    sol(5,v) = CrackSolutionGeneralCohesive_InsertVariables_Neri(CfVec(v)/1255,1e-8,Gammas(v),'L',3e-3,HandSadParab);
    cohsvMdlDef = @(x) (1 - x.^2).*(0.5*tanh(3*x + 1) + 0.5) + (x + 1).^2.*(0.5*tanh(-3*x - 1) + 0.5);
    sol(6,v) = CrackSolutionGeneralCohesive_InsertVariables_Neri(CfVec(v)/1255,1e-8,Gammas(v),'L',9e-3,cohsvMdlDef);

    
    %% plot
    CohesiveExmp = figure;
    for i=1:length(sol)
        subplot(length(CfVec),3,1+(v-1)*3);
        hold on;
        plot(Xchsv,sol(i,v).tauOfx(Xchsv));
        
        subplot(length(CfVec),3,2+(v-1)*3);
        hold on;
        plot(sol(i,v).x,sol(i,v).Ux);
        
        subplot(length(CfVec),3,3+(v-1)*3);
        hold on;
        plot(sol(i,v).x,sol(i,v).vx);
    end
    pause(0.5);
    drawnow
end