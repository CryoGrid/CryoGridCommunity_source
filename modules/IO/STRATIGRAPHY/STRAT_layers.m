classdef STRAT_layers
    
    properties
        depth
        variable_names
        variable_values
        variable_gridded
    end
    
    methods
        
        function xls_out = write_excel(strat)
            xls_out = {'STRATIGRAPHY','index',NaN,'both time-invariable parameters and initial state of time-variable parameters, add more cells if necessary',NaN,NaN,NaN,NaN,NaN;'STRAT_layers',1,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;'depth','waterIce','mineral','organic','Xice','field_capacity','soil_type','saltConc','deltaT';'[m]','[-]','[-]','[-]','[-]','[-]','code','mol/m3',NaN;'TOP',NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0.500000000000000,0.500000000000000,0,0,0.500000000000000,0,500,0;0.500000000000000,0.500000000000000,0.500000000000000,0,0,0.500000000000000,0,1000,1;10,0.0300000000000000,0.970000000000000,0,0,0.500000000000000,0,0,1;'BOTTOM',NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;'STRATIGRAPHY_END',NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN};
        end
        
        function strat = initalize_from_file(strat, section)
            pos_list = get_range_TOP_BOTTOM(section);
            strat.depth = cell2mat(section(pos_list(1,1):pos_list(1,2), 1));
            strat.variable_names={};
            strat.variable_values=[];
            i=2;
            field=cell2mat(section(pos_list(1,1)-3, i));
            
            while i<=size(section,2) && ~isnan(field(1))
                field=cell2mat(section(pos_list(1,1)-3, i));
                strat.variable_names=[strat.variable_names section{pos_list(1,1)-3, i}];
                strat.variable_values = [strat.variable_values cell2mat(section(pos_list(1,1):pos_list(1,2), i))];
                i=i+1;
            end
        end
        
        function strat = interpolate_to_grid(strat, grid)
            strat.variable_gridded = repmat(grid.MIDPOINTS .* 0, 1, size(strat.variable_values,2));
            for j=1:size(strat.variable_values,1)-1
                range = grid.MIDPOINTS > strat.depth(j,1) & grid.MIDPOINTS <= strat.depth(j+1,1);
                for i=1:size(strat.variable_values,2)
                    strat.variable_gridded(range,i) = strat.variable_gridded(range,i) + strat.variable_values(j,i);
                end
            end
            range = grid.MIDPOINTS > strat.depth(end,1);
            for i=1:size(strat.variable_values,2)
                    strat.variable_gridded(range,i) = strat.variable_gridded(range,i) + strat.variable_values(end,i);
            end
        end
        
    end
    
end
