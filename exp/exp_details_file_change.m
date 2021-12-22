function exp_details_file_change(date)
% date is should be cell array that contains strings with dates.
% example={'2011-07-13'}
%The function takes the exp_detail.txt from c:\frics and copies it to all experiments and the
%relevant directories at that date.

for k=1:length(date)
    
    exp_dir=dir(['c:\frics\' date{k}]);
    source='c:\frics\exp_details.txt';
    
    for j=1:length(exp_dir)
        
        if length(exp_dir(j).name)==8 && exp_dir(j).isdir==1
            destination=['c:\frics\' date{k} '\' exp_dir(j).name];
            copyfile(source,destination);
            
            destination=['c:\frics\' date{k} '\' exp_dir(j).name '\cal'];
            if exist(destination,'dir')==7
                copyfile(source,destination);
            end
            
        end
    end
end