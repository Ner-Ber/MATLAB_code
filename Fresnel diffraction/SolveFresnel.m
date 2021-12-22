%% define parameters
Lambda = 532*1e-9;
k = 2*pi/Lambda;
res = 97500;
dx = 0.05;
x_true = 0.5:dx:(frameLength+0.5);                % sub pixles coordinates
res_true = res/dx;
spatialVec_true = x_true./res;

scratchWidth = 0.0001;
S_W = [0.005 0.025 0.05 0.125 0.18];
stepWidth = S_W(1);
%% prepare initial conditions
Psi0 = @(p) ((tanh((mod(p/2/scratchWidth,1)-0.25)./stepWidth).*tanh((0.75-mod(p/2/scratchWidth,1))./stepWidth))+1).*0.5;
Psi0Damped = @(p) Psi0(p).*exp(-((p-mean(spatialVec_true)).^2)./(0.55e-6));


% spatialVec_true = -0.6e-3:8e-7:0.6e-3;
% Psi0 = @(p) ((tanh((p/2/scratchWidth+0.25)./stepWidth).*tanh((0.25-p/2/scratchWidth)./stepWidth))+1).*0.5;
% Nsteps = 5;
% D = (Nsteps-1)*scratchWidth;
% Dunit = 2*D/(Nsteps-1);
% if Nsteps==1
%     combLocation = 0;
% else
%     combLocation = -D:Dunit:D;
% end
% [~,combLocationPix] = min(abs(bsxfun(@minus,spatialVec_true(:),combLocation(:)')),[],1);
% comb = zeros(size(spatialVec_true));
% comb(combLocationPix) = 1;
% Psi0conv = @(p) conv(Psi0(p),comb,'same');
% Psi0Damped = @(p) Psi0conv(p);

%% solve integral
% Y = linspace(1e-4,0.008,5);
Y = 0.002;
xmin = 2*min(spatialVec_true)-max(spatialVec_true);
xmax = 2*max(spatialVec_true)-min(spatialVec_true);
I = nan(length(Y),length(spatialVec_true));
warning('off');
for j = 1:length(Y);
    yj = Y(j);
    for i = 1:length(spatialVec_true);
        x_i = spatialVec_true(i);
                Intgrnd = @(p) Psi0Damped(p).*sqrt(k./(2*pi*1i*yj)).*exp(1i.*k.*(yj+(x_i-p).^2./(2*yj)));
                I(j,i) = integral(Intgrnd,xmin,xmax);
%         Intgrnd = @(p) Psi0Damped(spatialVec_true).*sqrt(k./(2*pi*1i*yj)).*exp(1i.*k.*(yj+(x_i-spatialVec_true).^2./(2*yj)));
%         I(j,i) = trapz(spatialVec_true,Intgrnd(spatialVec_true));
    end
end
warning('on');

figure; hold on;
for i=1:length(Y);
    plot(spatialVec_true,abs(I(i,:)),'DisplayName',num2str(Y(i)))
end
plot(spatialVec_true,Psi0Damped(spatialVec_true),'DisplayName','original');