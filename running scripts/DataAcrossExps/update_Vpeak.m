load('Other\DAandREsidualsStruct.mat');
% load('DAandREsidualsStruct.mat');
D = dir_names
k1 = strfind(D,'.mat')
k2 = strfind(D,'allPhedi')
F = find(~cellfun(@isempty,k1) & ~cellfun(@isempty,k2))
maxSrchRegn = 0.03;
A_avgDist = 0.005;
dxIntrp = 5e-5;
intrpAx = -0.05:dxIntrp:0.05;
[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;


for j=1:length(F)
    %% iterate on events
    try
        disp(['loading ',D{F(j)}]);
        load(D{F(j)});
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
        relevantEvents = relevantEvents(relevantEvents~=1);
        
        peakVel_perPhd = cell(length(relevantEvents),1);
        peakVel_all_updated = nan(length(relevantEvents),1);
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
                relevantEvents = relevantEvents(relevantEvents~=1);
                
                
                PhediData = PhediStructCell{relevantEvents(iEv)}.PhediData;
                spaceAxis4Vel = PhediData.x_mins_x_tip_4vel;
                phediVel= PhediData.PhediVelocity;
                Cf = PhediData.Cf;
                %--- get peak velocity of mean phedi
                Sall = phedi_FWHMofMean(PhediStructCell{relevantEvents(iEv)},A_avgDist);
                VpeakXall = Sall.peakX;
                AvgPhediStructAll = Sall.AvgPhediStruct;
                [~,I] = min(abs(AvgPhediStructAll.X_mean- VpeakXall));
                peakVelFromMean = AvgPhediStructAll.Vel_mean_x(I);
                %--- get peak velocity of each phedi
%                 thisEvMaxVel = [];
%                 for p = 1:size(phediVel,2)
%                     nonInfNanLogic = ~(isinf(spaceAxis4Vel(:,p)) | isinf(phediVel(:,p)) | isnan(spaceAxis4Vel(:,p)) | isnan(phediVel(:,p)));
%                     intrpdVel = interp1(spaceAxis4Vel(nonInfNanLogic,p),phediVel(nonInfNanLogic,p),intrpAx);
%                     intrpdVelSmth = smooth(intrpdVel,round(A_avgDist/dxIntrp));
%                     maxVel = max(intrpdVelSmth(abs(intrpAx)<=maxSrchRegn));
%                     thisEvMaxVel = cat(1,thisEvMaxVel,maxVel);
%                 end
%                 peakVel_perPhd{iEv} = thisEvMaxVel;
                peakVel_all_updated(iEv) = peakVelFromMean;
                
            catch
                disp('');
            end
        end
        
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        J = arrayfun(@(S) strcmpi(S.Name,Name), myStruct);
        
%         myStruct(J).peakVel_perPhd = peakVel_perPhd;
        myStruct(J).peakVel_all_updated = peakVel_all_updated;
        
        
    catch
        disp([D{F(j)},' Failed to save']);
    end
end