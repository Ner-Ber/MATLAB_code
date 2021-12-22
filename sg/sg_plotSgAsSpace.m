function sg_plotSgAsSpace(PhediStruct)

SgData = PhediStruct.SgData;
% BigPicRotStruct = PhediStruct.BigPicRotStruct;
N = size(SgData.x_mins_x_tips,2);

%--- find t_tips at sg's
% x_sg = SgData.x_sg*1e-3; % sg location in [m]
% [~,minIdx] = min(abs(bsxfun(@minus,BigPicRotStruct.x(:),x_sg(:)')));
% t_tips = BigPicRotStruct.frontTime_interp(minIdx)/BigPicRotStruct.fps;
%--- create t-t_tips vectors
% t_mins_t_tips = bsxfun(@minus,SgData.t*1e-3,t_tips(:)');
%--- find the average velocity of front in some vicinity
% vicinity_m = 0.001;
% vicinity_pix = round(vicinity_m/BigPicRotStruct.res);
% [~,minIdx_vel] = min(abs(bsxfun(@minus,BigPicRotStruct.frontVelLoc_interpM(:),x_sg(:)')));
% ii = -vicinity_pix:vicinity_pix;
% pixRange = bsxfun(@plus,ii(:),minIdx_vel(:)');
% VelInVic = mean(BigPicRotStruct.frontVel_interpMperS(pixRange),1);
%--- create spacial vector:
% x_mins_x_tips = -bsxfun(@times,t_mins_t_tips,VelInVic(:)');


x_mins_x_tips = SgData.x_mins_x_tips;

figure; hold on;
MyCol = MyVaryColor(N);
% for i=1:N
%     plot(SgData.x_mins_x_tips(:,i),1e-3*(SgData.Uxx(:,i)-mean(SgData.Uxx(1:100,i))),...
%         '.-','Color',MyCol(i,:));
% end
for i=1:N
    plot(x_mins_x_tips(:,i),1e-3*(SgData.Uxx(:,i)-mean(SgData.Uxx(1:100,i))),...
        '.-','Color',MyCol(i,:));
end
xlabel('x-x_{tip} [m]'); ylabel('Uxx [strain]');
xlim([-0.05 0.05])
colormap(MyCol);
cb = colorbar;
cb.Ticks = movmean(linspace(0,1,N+1),2,'Endpoints','discard');
cb.TickLabels = mat2cell(round(SgData.x_sg*1e-3,3),1,N);
end