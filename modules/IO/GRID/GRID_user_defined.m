classdef GRID_user_defined

    properties
        GRID
        MIDPOINTS
        LAYERTHICK
        variable_names
        variable_gridded
    end
    
    methods
        
        function xls_out = write_excel(grid)
            xls_out = {'GRID','index',NaN;'GRID_user_defined',1,NaN;NaN,NaN,NaN;'upper','spacing','lower';'[m]','[m]','[m]';'TOP',NaN,NaN;0,0.0500000000000000,2;2,0.100000000000000,10;10,0.500000000000000,30;30,1,50;50,5,100;100,10,150;150,50,500;500,100,1000;'BOTTOM',NaN,NaN;'GRID_END',NaN,NaN};
        end
        
        function grid = initalize_from_file(grid, section)
            pos_list = get_range_TOP_BOTTOM(section);
            grid_breaks = cell2mat(section(pos_list(1,1):pos_list(1,2), 1:3));
            grid_breaks(2:end,1) = grid_breaks(2:end,1) + grid_breaks(2:end,2);
            grid.GRID =[];
            for i=1:size(grid_breaks)
                grid.GRID = [grid.GRID; [grid_breaks(i,1):grid_breaks(i,2):grid_breaks(i,3)]'];
            end
        end
        
        function grid = reduce_grid(grid, forcing)
            grid.GRID(grid.GRID > forcing.PARA.domain_depth)=[]; 
            grid.MIDPOINTS = (grid.GRID(2:end,1) + grid.GRID(1:end-1,1))./2;
            grid.LAYERTHICK = (grid.GRID(2:end,1) - grid.GRID(1:end-1,1));
        end
    end
    
end

