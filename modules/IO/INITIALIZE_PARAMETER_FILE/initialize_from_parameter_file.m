function [class_list, FIRST, LAST, class_strat] = initialize_from_parameter_file(parameter_file)
    %addpath('modules')
    %addpath('modules/initialize_from_parameter_file')
    addpath('../')
    [dummy1, dummy2, status_info] = xlsread(parameter_file);
    
    CONST_file = 'CONSTANTS_excel.xlsx';
    [dummy1, dummy2, CONST_info] = xlsread(CONST_file);
    
    class_list={};
    FIRST=1;
    LAST=[];
    

    
    [lat, lon, alt, domain_depth] = read_coordinates(status_info);
    grid = get_grid(status_info);
    strat = get_stratigraphy(status_info);
    T_initial = get_T_initial(status_info);
    class_strat = get_class_strat(status_info);
    
    

    class_positions = get_range(status_info, 'CLASS');

    for i=1:size(class_positions,1)
        %create and initialize all classes
        class_handle = str2func(status_info{class_positions(i,1)+1,1});
        class_index = status_info{class_positions(i,1)+1,2};
        class=class_handle();
        class=initialize(class); %set default values for each class
        
        %read and overwrite CONST for each class
        class_variables = fieldnames(class.CONST);
        for j=1:size(class_variables,1)
            for k=1:size(CONST_info,1)
               if strcmp(CONST_info{k,1}, class_variables(j,1))
                  class.CONST.(class_variables{j}) = CONST_info{k,2};
               end
            end
        end
        
        %read and overwrite default PARA for all classes
        class_variables = fieldnames(class.PARA);
        search_window = status_info(class_positions(i,1)+1:class_positions(i,2)-1,:);
        
        for j=1:size(class_variables,1)
            for k=1:size(search_window,1)
               if strcmp(search_window{k,1}, class_variables(j,1))
                  class.PARA.(class_variables{j}) = search_window{k,2};
               end
            end
        end
        class_list=[class_list; {class class_index}];
        
        
        %initialize STATVAR and additional parameters and produce the first
%         %stratigraphy
%        1. grid
%        2. altitude
%        3. T_initial
%        4. stratigraphy
%         
        
    end
    
end



function [lat, lon, alt, domain_depth] = read_coordinates(status_info)

pos_list = get_range(status_info, 'COORDINATES');
section = status_info(pos_list(1,1):pos_list(1,2),:);

for i=1:size(section,1)
    if strcmp(section{i,1}, 'latitude')
        lat=section{i,2};
    end
    if strcmp(section{i,1}, 'longitude')
        lon=section{i,2};
    end
    if strcmp(section{i,1}, 'surface_altitude')
        alt=section{i,2};
    end
    if strcmp(section{i,1}, 'domain_depth')
        domain_depth=section{i,2};
    end
end

end


function grid = get_grid(status_info)

pos_list = get_range(status_info, 'GRID');
section = status_info(pos_list(1,1):pos_list(1,2),:);
pos_list = get_range_TOP_BOTTOM(section);

grid_breaks = cell2mat(section(pos_list(1,1):pos_list(1,2), 1:3));
grid_breaks(2:end,1) = grid_breaks(2:end,1) + grid_breaks(2:end,2);

grid =[];
for i=1:size(grid_breaks)
    grid = [grid; [grid_breaks(i,1):grid_breaks(i,2):grid_breaks(i,3)]'];
end

end


function strat = get_stratigraphy(status_info)

pos_list = get_range(status_info, 'STRATIGRAPHY');
section = status_info(pos_list(1,1):pos_list(1,2),:);
pos_list = get_range_TOP_BOTTOM(section);

name_pos = pos_list(1,1)-3;

strat{1}.name = 'depth';
strat{1}.value = cell2mat(section(pos_list(1,1):pos_list(1,2), 1));

j=2; 
while ~isnan(section{name_pos,j})
    strat{j}.name = section{name_pos,j};
    strat{j}.value = cell2mat(section(pos_list(1,1):pos_list(1,2), j));
    j=j+1;
end
end


function T_initial = get_T_initial(status_info)

pos_list = get_range(status_info, 'T_INITIAL');
section = status_info(pos_list(1,1):pos_list(1,2),:);
pos_list = get_range_TOP_BOTTOM(section);

T_initial = cell2mat(section(pos_list(1,1):pos_list(1,2), 1:2));
end


function class_strat = get_class_strat(status_info)

pos_list = get_range(status_info, 'CLASS_STRATIGRAPHY');
section = status_info(pos_list(1,1):pos_list(1,2),:);
pos_list = get_range_TOP_BOTTOM(section);

class_strat.upper_depths = cell2mat(section(pos_list(1,1):pos_list(1,2), 1));
class_strat.class_names = section(pos_list(1,1):pos_list(1,2), 2);
class_strat.index = cell2mat(section(pos_list(1,1):pos_list(1,2), 3));
end


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

end


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

end