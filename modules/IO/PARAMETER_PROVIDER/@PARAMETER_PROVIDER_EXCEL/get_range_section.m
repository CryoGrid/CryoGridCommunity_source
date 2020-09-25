function pos_list = get_range_section(self, section)
% GET_RANGE_SECTION Get range of cells in cell array that make up the section
%   and the individual classes defined in the section.

%   ARGUMENTS:
%   section: a string specifying which section to list (e.g. 'FORCING' or 'OUT')
%
%   RETURNS: 
%   pos_list: and array containing the positions of the first and last
%   lines of the section

pos_list=[];

for i=1:size(self.config_data,1)
    if strcmp(self.config_data{i,1}, section)
        j=i+1;
        while ~contains(num2str(self.config_data{j,1}), '_END')
            % disp(self.config_data{j,1});
            j=j+1;
        end
        pos_list=[pos_list; i j];
    end
end
end

