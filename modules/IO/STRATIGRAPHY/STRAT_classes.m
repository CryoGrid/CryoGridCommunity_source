classdef STRAT_classes

    
    properties
        depth
        class_name
        class_index
        snow_class
    end
    
    methods
        function strat = initalize_from_file(strat, section)
            for i=1:size(section,1)
                if strcmp(section{i,1}, 'snow_class')
                    strat.snow_class.name = section{i,2};
                    strat.snow_class.index = section{i,3};
                end
            end
            pos_list = get_range_TOP_BOTTOM(section);
            strat.depth = cell2mat(section(pos_list(1,1):pos_list(1,2), 1));
            strat.class_name =  section(pos_list(1,1):pos_list(1,2), 2);
            strat.class_index = cell2mat(section(pos_list(1,1):pos_list(1,2), 3));
        end
        
        function strat = interpolate_to_grid(strat, grid)
            depth2=strat.depth;
            for i=1:size(depth2)-1
                depth2(i,1) = strat.depth(i+1,1);
            end
            depth2(end,1) = grid.GRID(end,1);
            strat.depth = [strat.depth depth2];
        end
        
    end
    
end

