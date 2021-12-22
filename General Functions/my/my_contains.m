function containLogical = my_contains(C,STR)

containLogical = ~cellfun(@isempty,strfind(C, STR));

end