function my_legend_add(legendV)
% The function adds legends to figure which already has some legends.
%legendV should be char vector column oriented or cell array of characters


axe=findobj(gca,'Type','line');
if (ischar(legendV))
legendV=cellstr(legendV);
end
k=1;
    for j=length(axe):-1:1
        if   (strcmp(get(axe(j),'DisplayName'),''))
            set(axe(j),'DisplayName',legendV{k});
            k=k+1;
            if(k>length(legendV))
                break;
            end
        end
    end




