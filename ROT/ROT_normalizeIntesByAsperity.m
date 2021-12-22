function RowOverTime_normalized = ROT_normalizeIntesByAsperity(RowOverTime,MarkedMatCell)
    % RowOverTime_normalized = ROT_normalizeIntesByAsperity(RowOverTime,MarkedMatCell)
    %
    % RowOverTime - regular time space used in many functions
    % MarkedMatCell - a product of the function 'Movie_followAsperities'. A
    % cell array containing logical arrays the size of 'RowOverTime' where the
    % true imply on the asperity, and falsethe rest. Cell size is 1xnum. of
    % asperities.
    
    MarkedMatCellLogical = cellfun(@logical,MarkedMatCell,'UniformOutput',0);
    allLogical = MarkedMatCellLogical{1};
    RowOverTime_normalized_t = RowOverTime;
    for i=1:length(MarkedMatCellLogical)
        AsperRowInt = trapz(RowOverTime.*MarkedMatCellLogical{i},2);
        thisAsperNormlzd = bsxfun(@rdivide,RowOverTime.*MarkedMatCellLogical{i},AsperRowInt);
        RowOverTime_normalized_t(MarkedMatCellLogical{i}) = thisAsperNormlzd(MarkedMatCellLogical{i});
        allLogical = allLogical | MarkedMatCellLogical{i};
    end
    
    %-- eliminate all pixels outside asperities
    RowOverTime_normalized = RowOverTime_normalized_t;
    RowOverTime_normalized(~allLogical) = 0;
    
end