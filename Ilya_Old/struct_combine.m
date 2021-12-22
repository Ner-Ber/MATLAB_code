function struct=struct_combine(strct1,strct2)

fName=fieldnames(strct1);

for j=1:length(fName)
    struct.(fName{j})=[strct1.(fName{j}) strct2.(fName{j})]
end