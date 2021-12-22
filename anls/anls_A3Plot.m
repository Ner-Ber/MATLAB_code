function anls_A3Plot(anl,kind)
% unlike anls_A this function doesn't take average values from different
% strain gages.

smtX=5;
x_Start=80;
x_End=130;

%figure;
hold all;
for j=1:length(anl.event)
    
    if strcmp(kind,'Rayleigh')
        %SubRayleigh
        anl.front{j}.v(anl.front{j}.v>1300)=NaN;
    else
        %SuperShear
        anl.front{j}.v=smooth(anl.front{j}.v,5);
        anl.front{j}.v(anl.front{j}.v<1300)=NaN;
    end
    
    Xc=my_smooth2(anl.front{j}.x1,anl.front{j}.x1-anl.front{j}.x2,smtX);
    Xc(Xc>10)=NaN;
        
    [~,index_Start]=min(abs(anl.front{j}.x1-x_Start));
    [~,index_End]=min(abs(anl.front{j}.x1-x_End));
    
    plot(anl.front{j}.v(index_Start:index_End),Xc(index_Start:index_End),'.-');
    my_legend_add(anl.lgnd{j});
    title([num2str(x_Start) '<x<' num2str(x_End)]);
    %ylim([0 5]);
    
end


end
