%% compare with Ilya's Xc
% relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
% relevantEvents = relevantEvents(relevantEvents~=1);
relevantEvents=18;
for iEv = relevantEvents
    try
        %% My Method
        findLocs = linspace(0.01,0.14,200);
        ChsvFromCntct_mega = struct([]);
        for f=1:length(findLocs)
%             ChsvFromCntct = ROT_findXcFromContact_fromFrame(PhediStructCell{iEv}.BigPicRotStruct,findLocs(f));
%             ChsvFromCntct = ROT_findXcFromContact_fromPix(PhediStructCell{iEv}.BigPicRotStruct,findLocs(f));
            ChsvFromCntct = ROT_findXcFromContact_fromPix(BigPicRotStruct,findLocs(f));
            ChsvFromCntct_mega = cat(1,ChsvFromCntct_mega,ChsvFromCntct);
        end
        Xloc = arrayfun(@(S) S.loc,ChsvFromCntct_mega);
        Xc_vec = arrayfun(@(S) S.Xc,ChsvFromCntct_mega);
        
        
        %% Ilya
        
        phE_Myn = struct;
        phE_Myn.Date = PhediStructCell{iEv}.ExperimentData.ExpDate;
        phE_Myn.Exp = PhediStructCell{iEv}.ExperimentData.ExpHour;
        phE_Myn.EventNum = PhediStructCell{iEv}.ExperimentData.eventNum;
        phE_Myn.t_trigger = [];
        phE_Myn.t= BigPicRotStruct.t;
        phE_Myn.lines= BigPicRotStruct.DataMat;
        phE_Myn.Included= ones(4,1280);
        phE_Myn.x= BigPicRotStruct.x;
        phE_Myn.firstLine= BigPicRotStruct.DataMat(1,:);
        phE_Myn.t= BigPicRotStruct.t*1e3;
        phE_Myn.x= BigPicRotStruct.x*1e3;
        
        
        [frontRaw front]=calc_frontXT_from_XAN(phE_Myn,-3,3,7,10,140);
        [frontRaw.xv frontRaw.v]=calc_front_v_standalone(frontRaw.x,frontRaw.t,5);
        [front.xv front.v]=calc_front_v_standalone(front.x1,front.t,5);
        
        x=front.x1; %another option is by front.x1
        t=front.t;
        for j=1:length(t)
            [~, index_vec(j)]=min(abs(phE_Myn.t-t(j)));
        end
        phECut.t=phE_Myn.t(index_vec); % this gives a  vector of t-t_tip
        phECut.lines=phE_Myn.lines(index_vec,:)./repmat(phE_Myn.firstLine,length(phECut.t),1);
        phECut.frontX=x';
        phECut.x=phE_Myn.x;
        phECut.xOffset=repmat(phE_Myn.x,length(phECut.t),1)-repmat(phECut.frontX,1,length(phECut.x));
        dx=5;
        phECut.dx=dx;
        phECut.xf=-25;
        [~,index1]=min(abs(phECut.xOffset-(phECut.xf-dx)));
        [~,index2]=min(abs(phECut.xOffset-(phECut.xf+dx)));
        for j=1:length(index1)
            phECut.Af_from_t(j)=mean(phECut.lines(index2(j):index1(j),j),1);
        end
        phECut.xf2=-15;
        [~,index1]=min(abs(phECut.xOffset-(phECut.xf2-dx)));
        [~,index2]=min(abs(phECut.xOffset-(phECut.xf2+dx)));
        for j=1:length(index1)
            phECut.Af2_from_t(j)=mean(phECut.lines(index2(j):index1(j),j),1);
        end
        phECut.xf3=-35;
        [~,index1]=min(abs(phECut.xOffset-(phECut.xf3-dx)));
        [~,index2]=min(abs(phECut.xOffset-(phECut.xf3+dx)));
        for j=1:length(index1)
            phECut.Af3_from_t(j)=mean(phECut.lines(index2(j):index1(j),j),1);
        end
        clear index1 index2;
        %---calc Xc
        lines=phECut.lines-repmat(phECut.Af2_from_t,size(phECut.xOffset,1),1);
        lines=lines./(1-repmat(phECut.Af2_from_t,size(phECut.xOffset,1),1));
        [~,index1]= min(abs(lines-0.37));
        for j=1:size(phECut.xOffset,2)
            phECut.Xc(j)=phECut.xOffset(index1(j),j);
        end
        clear index1;
        %------Find Af from space .
        [~,index1]=min(abs(phECut.xOffset'-(phECut.xf2-dx)));
        [~,index2]=min(abs(phECut.xOffset'-(phECut.xf2+dx)));
        for j=1:length(index1)
            phECut.Af2_from_x(j)=mean(phECut.lines(j,index1(j):index2(j)),2);
        end
        
        %% new myn
        chsvFromContact = ROT_findXcFromContact(BigPicRotStruct);
        
        %% plot both
        figure; IDT_PlotRowOverTime(BigPicRotStruct);
        hold on;
        plot(frontRaw.x*1e-3,frontRaw.t*1e-3,'g','LineWidth',1.5);
        
        
        figure; hold on;
        plot(phECut.x,abs(phECut.Xc),'.-','DisplayName','Ilya');
        hold on;
        plot(Xloc*1e3,abs(Xc_vec*1e3),'.-','DisplayName','Old');
        plot(BigPicRotStruct.x*1e3,abs(chsvFromContact.Xc*1e3),'.-','DisplayName','New');
    catch
    end
end