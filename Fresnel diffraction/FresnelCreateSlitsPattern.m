function [InterferenceProfile, Psi0Damped] = FresnelCreateSlitsPattern(spatialVec, varargin)
%[InterferenceProfile, Psi0Damped] = FresnelCreateSlitsPattern(spatialVec,
%Y, scratchWidth, stepWidth, Lambda, sigmaDamp)

%% define parameters
[Y, scratchWidth, stepWidth, Lambda, sigmaDamp] = setDefaults4function(varargin,...
    0.002, 1e-4, 0.005, 532*1e-9, sqrt(0.55e-6./2));

k = 2*pi/Lambda;



%% prepare initial conditions
Psi0 = @(p) ((tanh((mod(p/2/scratchWidth,1)-0.25)./stepWidth).*tanh((0.75-mod(p/2/scratchWidth,1))./stepWidth))+1).*0.5;
% Psi0Damped = @(p) Psi0(p).*exp(-((p-mean(spatialVec)).^2)./(2*sigmaDamp^2));
Psi0Damped = @(p) Psi0(p);


% spatialVec = -0.6e-3:8e-7:0.6e-3;
% Psi0 = @(p) ((tanh((p/2/scratchWidth+0.25)./stepWidth).*tanh((0.25-p/2/scratchWidth)./stepWidth))+1).*0.5;
% Nsteps = 5;
% D = (Nsteps-1)*scratchWidth;
% Dunit = 2*D/(Nsteps-1);
% if Nsteps==1
%     combLocation = 0;
% else
%     combLocation = -D:Dunit:D;
% end
% [~,combLocationPix] = min(abs(bsxfun(@minus,spatialVec(:),combLocation(:)')),[],1);
% comb = zeros(size(spatialVec));
% comb(combLocationPix) = 1;
% Psi0conv = @(p) conv(Psi0(p),comb,'same');
% Psi0Damped = @(p) Psi0conv(p);

%% solve integral
xmin = 2*min(spatialVec)-max(spatialVec);
xmax = 2*max(spatialVec)-min(spatialVec);
InterferenceProfile = nan(length(Y),length(spatialVec));
warning('off');
for j = 1:length(Y);
    for i = 1:length(spatialVec);
        x_i = spatialVec(i);
                Intgrnd = @(p) Psi0Damped(p).*sqrt(k./(2*pi*1i*Y)).*exp(1i.*k.*(Y+(x_i-p).^2./(2*Y)));
                InterferenceProfile(j,i) = integral(Intgrnd,xmin,xmax);
%         Intgrnd = @(p) Psi0Damped(spatialVec).*sqrt(k./(2*pi*1i*yj)).*exp(1i.*k.*(yj+(x_i-spatialVec).^2./(2*yj)));
%         I(j,i) = trapz(spatialVec,Intgrnd(spatialVec));
    end
end
warning('on');


end