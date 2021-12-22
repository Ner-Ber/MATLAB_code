function DataStructUpdated = phedi_add_velocity_to_struct(DataStruct)
% DataStructUpdated = phedi_add_velocity_to_struct(DataStruct)
%
% phedi_add_velocity_to_struct will add 

PhediLocation = DataStruct.PhediData.PhediLocation;
timeVec = DataStruct.PhediData.timeVec;
t_mins_t_tip = DataStruct.PhediData.t_mins_t_tip;
x_mins_x_tip = DataStruct.PhediData.x_mins_x_tip;

[PhediVelocity, timeVec_4vel] = phedi_calcVelocity(PhediLocation,timeVec);
t_mins_t_tip_4vel = conv2(t_mins_t_tip,[0.5 0.5]','valid');
x_mins_x_tip_4vel = conv2(x_mins_x_tip,[0.5 0.5]','valid');

DataStructUpdated = DataStruct;
DataStructUpdated.PhediData.PhediVelocity = PhediVelocity;
DataStructUpdated.PhediData.timeVec_4vel = timeVec_4vel;
DataStructUpdated.PhediData.t_mins_t_tip_4vel = t_mins_t_tip_4vel;
DataStructUpdated.PhediData.x_mins_x_tip_4vel = x_mins_x_tip_4vel;

end