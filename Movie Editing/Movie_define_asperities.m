function [vallyPairsCell,asper_num]=Movie_define_asperities(RowOverTimeCell,Rows,factor,min_gap,asperity_min_size)

asper_num=zeros(1,length(Rows));

vallyPairsCell={};
for rowIdx = 1:length(Rows)
    
    First_line=RowOverTimeCell{rowIdx}(1,:);
%         figure; hold all
%         plot(First_line)
    
    
    threshold=factor*mean(First_line);
%         line([0 2000],[threshold threshold],'Color','k');
    logic=(First_line>threshold);
    a=find(logic);
    
    
    k=1;
    V=[];
    V(1,1)=a(1);
    for i=1:length(a)-1
        if a(i)<a(i+1)-min_gap
            V(k,2)=a(i);
            V(k+1,1)=a(i+1);
            k=k+1;
        end
    end
    V(k,2)=a(end);
    
    
    Val=[];
    for i=1:size(V,1)
        if V(i,2)-V(i,1) > asperity_min_size
            Val=cat(1,Val,V(i,:));
        end
    end
    asper_num(rowIdx)=size(Val,1);
    vallyPairsCell{rowIdx}=Val;
end

end