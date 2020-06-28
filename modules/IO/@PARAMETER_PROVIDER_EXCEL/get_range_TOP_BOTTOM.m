function pos_list = get_range_TOP_BOTTOM(section_data)

pos_list=[];

for i=1:size(section_data,1)
    if strcmp(section_data{i,1}, 'TOP')
        j=i+1;
        while ~strcmp(section_data{j,1}, 'BOTTOM')
            j=j+1;
        end
        pos_list=[i+1 j-1];
    end
end
