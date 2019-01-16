function pos_list = get_range_TOP_BOTTOM(status_info)

pos_list=[];

for i=1:size(status_info,1)
    if strcmp(status_info{i,1}, 'TOP')
        j=i+1;
        while ~strcmp(status_info{j,1}, 'BOTTOM')
            j=j+1;
        end
        pos_list=[i+1 j-1];
    end
end
