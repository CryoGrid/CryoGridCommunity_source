function pos_list = get_range(self, keyword)

pos_list=[];

for i=1:size(self.config_data,1)
    if strcmp(self.config_data{i,1}, keyword)
        j=i+1;
        while ~strcmp(self.config_data{j,1}, [keyword '_END'])
            % disp(self.config_data{j,1});
            j=j+1;
        end
        pos_list=[pos_list; i j];
    end
end
end
