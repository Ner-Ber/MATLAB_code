function PhediStructUpdated = phedi_smoothWsgolay(DataStruct,varargin)
    % PhediStructUpdated = phedi_smoothWsgolay(PhediStruct,order, frameLength)
    %
    %       OPTIONAL (DEFAULT)
    %       order(6)
    %       frameLength(29)
    %
    %
    %   PhediStructUpdated will use Savitzky-Golay smoothing and
    %   differentiation to update the location and velocity in the main
    %   structure.
    %   The output is a similar structure with the old field saved under a
    %   new name and a new field instead
    
    %% set defaults and save old data
    [order, framelen, derivOrdr] = setDefaults4function(varargin,6,25,1);
    
    PhediStructUpdated = DataStruct;
    %-- if there's already a smoothed version:
    if isfield(DataStruct,'PhediDataPreSmth')
        PhediStructUpdated.PhediDataSimpSmth = PhediStructUpdated.PhediData;
        PreSmthStruct = PhediStructUpdated.PhediDataPreSmth;
    else
        PreSmthStruct = PhediStructUpdated.PhediData;
        PhediStructUpdated.PhediDataPreSmth = PreSmthStruct;
    end
    
    %% smooth phedis
    [PhediLocation, timeVec] = deal(PreSmthStruct.PhediLocation,PreSmthStruct.timeVec);
    [b, g] = sgolay(order,framelen);
    cntrKernel = b((framelen+1)/2,:);
%     PhediLocation_center = conv(PhediLocation2,b((framelen+1)/2,:),'valid');
    PhediLocation_center = conv2(PhediLocation,cntrKernel(:),'valid');
    
    
    bgnKernel = b(end:-1:(framelen+3)/2,:);
    PhediLocation_begin = bgnKernel*PhediLocation(framelen:-1:1,:);
    endKernel = b((framelen-1)/2:-1:1,:);
    PhediLocation_end = endKernel*PhediLocation(end:-1:end-(framelen-1),:);    
    PhediLocation_smth = [PhediLocation_begin; PhediLocation_center; PhediLocation_end];
    
    %% compute velocity
    derivKrnl = factorial(derivOrdr)/(-(timeVec(2)-timeVec(1)))^derivOrdr * g(:,derivOrdr+1);
    PhediVelocity = conv2(PhediLocation, derivKrnl(:), 'same');
    
    %% add data to structure
    SgSmthStruct = PreSmthStruct;
    SgSmthStruct.PhediLocation = PhediLocation_smth;
    SgSmthStruct.PhediVelocity = PhediVelocity;
    SgSmthStruct.timeVec_4vel = timeVec;
    SgSmthStruct.t_mins_t_tip_4vel = SgSmthStruct.t_mins_t_tip;
    SgSmthStruct.x_mins_x_tip_4vel = SgSmthStruct.x_mins_x_tip;
    
    PhediStructUpdated.PhediDataSGsmth = SgSmthStruct;
    
end