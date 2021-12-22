function prof_parameters = profilometer_readParameters(filename)

%% open file to read
fid1=fopen(filename,'r');

%% organize data
%---Read Header Information
parametersReadOut=textscan(fid1,'%s','delimiter',{':','='});
parametersNames = {'Pen','Freq','distanceX','StepDistX','VelocityX','NumOfRows','RowStep','VelocityY'};

parametersVals = parametersReadOut{1}(2:2:end);
parametersNumVals = cellfun(@(x) sscanf(x,'%f'),parametersVals,'UniformOutput',0);
A = cellfun(@(x) regexp(x,'\D*','Match'),parametersVals,'UniformOutput',0);
A{6} = {nan};
parametersUnits = cellfun(@(x) x{end},A,'UniformOutput',0);

%--- transform to struct
fieldNames = {'Val','unit'};
ParametersData = cell2struct([parametersNumVals,parametersUnits],fieldNames,2);
prof_parameters = struct(parametersNames{1},ParametersData(1),...
    parametersNames{2},ParametersData(2),...
    parametersNames{3},ParametersData(3),...
    parametersNames{4},ParametersData(4),...
    parametersNames{5},ParametersData(5),...
    parametersNames{6},ParametersData(6),...
    parametersNames{7},ParametersData(7),...
    parametersNames{8},ParametersData(8));


%% done 
fclose(fid1);



end