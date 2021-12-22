function PhediStructUpdated = phedi_smoothByfreqFilter(DataStruct,varargin)
    % PhediStructUpdated = phedi_smoothWsgolay(PhediStruct,order, frameLength)
    %
    %       OPTIONAL (DEFAULT)
    %       cutOff(fps/5)
    %
    %
    %   PhediStructUpdated will use Savitzky-Golay smoothing and
    %   differentiation to update the location and velocity in the main
    %   structure.
    %   The output is a similar structure with the old field saved under a
    %   new name and a new field instead
    
    %% set defaults
    Fs = DataStruct.CamMeta.FrameRate;
    [cutOff] = setDefaults4function(varargin,Fs/5);
    
    
    PhediStructUpdated = DataStruct;
    %-- if there's already a smoothed version:
    if isfield(DataStruct,'PhediDataPreSmth')
        PhediStructUpdated.PhediDataSimpSmth = PhediStructUpdated.PhediData;
        PreSmthStruct = PhediStructUpdated.PhediDataPreSmth;
    else
        PreSmthStruct = PhediStructUpdated.PhediData;
        PhediStructUpdated.PhediDataPreSmth = PreSmthStruct;
    end
    %%
%     time = PreSmthStruct.timeVec;
%     L = length(time);
%     LocationSpectrm = fft(PreSmthStruct.PhediLocation);
%     LocationSpectrmPix = fft(PreSmthStruct.PhediLocationPix);
%     LocationSpectrm_shifted = fftshift(LocationSpectrm);
%     LocationSpectrmPix_shifted = fftshift(LocationSpectrmPix);
%     f = (-L/2:L/2-1)*(Fs/L);
%     LocationSpectrm_zero = LocationSpectrm_shifted;
%     LocationSpectrmPix_zero = LocationSpectrmPix_shifted;
%     
%     LocationSpectrm_zero(abs(f)>=cutOff,:) = 0;
%     LocationSpectrmPix_zero(abs(f)>=cutOff,:) = 0;
%     
%     PhediLocation_filtered = ifft(ifftshift(LocationSpectrm_zero));
%     PhediLocationPix_filtered = ifft(ifftshift(LocationSpectrmPix_zero));
    
    %%
    kL = round(Fs/cutOff);
    kernel = ones(kL,1)./kL;
    PhediLocation_filtered = conv2(PreSmthStruct.PhediLocation,kernel,'same');
    PhediLocationPix_filtered = conv2(PreSmthStruct.PhediLocationPix,kernel,'same');
    
    %% update new structure
    PhediStructUpdated.PhediData.PhediLocation = PhediLocation_filtered;
    PhediStructUpdated.PhediData.PhediLocationPix = PhediLocationPix_filtered;
    
    
    %% add new phedi velocity
    PhediStructUpdated = phedi_add_velocity_to_struct(PhediStructUpdated);
end