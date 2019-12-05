classdef STRAT_linear
    
    properties
        depth
        variable_names
        variable_values
        variable_gridded
    end
    
    methods
        
        function xls_out = write_excel(strat)
            xls_out = {'STRATIGRAPHY','index';'STRAT_linear',1;NaN,NaN;'depth','T';'[m]','[degree C]';'TOP',NaN;0,1;1,0;10,-5;100,0;5000,20;'BOTTOM',NaN;'STRATIGRAPHY_END',NaN};
        end
        
        function strat = initalize_from_file(strat, section)
            pos_list = get_range_TOP_BOTTOM(section);
            strat.depth = cell2mat(section(pos_list(1,1):pos_list(1,2), 1));
            strat.variable_names={};
            strat.variable_values=[];
            i=2;
            field=cell2mat(section(pos_list(1,1)-3, i));
            while i<=size(section,2) && ~isnan(field(1))
                strat.variable_names=[strat.variable_names section{pos_list(1,1)-3, i}];
                strat.variable_values = [strat.variable_values cell2mat(section(pos_list(1,1):pos_list(1,2), i))];
                i=i+1;
                field=cell2mat(section(pos_list(1,1)-3, i));
            end
        end
        
       function strat = interpolate_to_grid(strat, grid)
           strat.variable_gridded = [];
           size(strat.variable_values,2);
           
            for i=1:size(strat.variable_values,2)
                    strat.variable_gridded = [strat.variable_gridded; interp1(strat.depth, strat.variable_values(:,i), grid.MIDPOINTS, 'linear')];
            end
        end
        
    end
    
end

