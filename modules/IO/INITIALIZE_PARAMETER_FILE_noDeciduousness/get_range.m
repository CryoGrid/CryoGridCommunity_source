function pos_list = get_range(status_info, keyword)

pos_list=[];

for i=1:size(status_info,1)
    if strcmp(status_info{i,1}, keyword)
        j=i+1;
        while ~strcmp(status_info{j,1}, [keyword '_END'])
            j=j+1;
        end
        pos_list=[pos_list; i j];
    end
end

