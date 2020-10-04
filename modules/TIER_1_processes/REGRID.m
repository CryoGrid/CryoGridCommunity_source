%========================================================================
% CryoGrid TIER1 library class for functions related to regridding operations
% S. Westermann, October 2020
%========================================================================

classdef REGRID < BASE
    
    
    methods
        function ground = regrid_split_merge(ground, extensive_variables) %ATTENTION: does not work stable!!
            %possibly better use regrid_split for short timesteps and once
            %in a while do regrid_full
            top_pos = ground.STATVAR.top_depth_rel2groundSurface;
            target_grid = ground.PARA.target_grid;
            target_layerThick = ground.PARA.target_layerThick;
            %first cell
            i=1;
            top_pos = top_pos + ground.STATVAR.layerThick(i,1);
            % pos_in_grid = max(1, sum(double(top_pos>target_grid))-1);  %grid starts with 0m
            pos_in_grid = max(1,min(sum(double(top_pos>target_grid)), size(target_layerThick,1)));
            if ground.STATVAR.layerThick(i,1) >= 1.5 .* target_layerThick(pos_in_grid,1)
                split_fraction = (ground.STATVAR.layerThick(i,1) - target_layerThick(pos_in_grid,1)) ./ ground.STATVAR.layerThick(i,1);
                ground = split_cell_extensive(ground, i, split_fraction, extensive_variables);
                top_pos = top_pos - ground.STATVAR.layerThick(i+1,1);
            elseif ground.STATVAR.layerThick(i,1) < 0.5 .* target_layerThick(pos_in_grid,1) %merge with previous cell
                top_pos = top_pos - ground.STATVAR.layerThick(i,1);
                ground = merge_cells_extensive(ground, i, i+1, extensive_variables);
                top_pos = top_pos + ground.STATVAR.layerThick(i,1);
            end
            
            %middle cells
            i=i+1;
            while i <= size(ground.STATVAR.layerThick,1)
                top_pos = top_pos + ground.STATVAR.layerThick(i,1);
                %pos_in_grid = max(1,sum(double(top_pos>target_grid))-1);  %grid starts with 0m
                pos_in_grid = max(1,min(sum(double(top_pos>target_grid)), size(target_layerThick,1)));
                if ground.STATVAR.layerThick(i,1) >= 1.5 .* target_layerThick(pos_in_grid,1)
                    %split_fraction = (ground.STATVAR.layerThick(i,1) - target_layerThick(pos_in_grid,1)) ./ ground.STATVAR.layerThick(i,1);
                    split_fraction =0.5;
                    ground = split_cell_extensive(ground, i, split_fraction, extensive_variables);
                    top_pos = top_pos - ground.STATVAR.layerThick(i+1,1);
                    %i=i+1;
                elseif ground.STATVAR.layerThick(i,1) < 0.5 .* target_layerThick(pos_in_grid,1) %merge with next cell
                    top_pos = top_pos - ground.STATVAR.layerThick(i,1);
                    last_cell = double(i == size(ground.STATVAR.layerThick,1));
                    ground = merge_cells_extensive(ground, i-last_cell, i+1-last_cell, extensive_variables);
                    top_pos = top_pos + ground.STATVAR.layerThick(i-last_cell,1);
                    
                elseif ground.STATVAR.layerThick(i-1,1)+ground.STATVAR.layerThick(i,1) < 1.5 .* target_layerThick(pos_in_grid,1)
                    
                    ground = merge_cells_extensive(ground, i-1, i, extensive_variables);
                    i=i-1;
                end
                i=i+1;
            end
            %the case that only one cell is left that must be annihilated must
            %be handled before, i.e. the class must be destroyed, whatever
            if ~strcmp(class(ground.NEXT), 'Bottom')  %set top_pos of next class to new depth
                ground.NEXT.STATVAR.top_depth_rel2groundSurface = top_pos;
            end
        end
        
        %split only
        function ground = regrid_split(ground, variable_list) %simple split function, merges cells if minimum thickness is reached
            top_pos = ground.STATVAR.top_depth_rel2groundSurface;
            min_thickness = ground.PARA.target_layerThick(1,1);

            regrid = find(ground.STATVAR.layerThick(:,1) < 0.5 .* min_thickness);
            
            while ~isempty(regrid)
                first_cell = double(regrid(end,1)==1);
                ground = merge_cells_extensive(ground, regrid(end,1)-1+first_cell, regrid(end,1)+first_cell, variable_list);
                regrid(end,:) =[];
            end
            
            %the case that only one cell is left that must be annihilated must
            %be handled before, e.g. the class must be removed
            if ~strcmp(class(ground.NEXT), 'Bottom')  %set top position of next class to new depth
                ground.NEXT.STATVAR.top_depth_rel2groundSurface = ground.STATVAR.top_depth_rel2groundSurface + sum(ground.STATVAR.layerThick,1);
            end
        end
        
        %complete regrid
        function ground = regrid_full(ground, variable_list) %regrids to the orginal grid  provided in initialization
            top_pos = ground.STATVAR.top_depth_rel2groundSurface;
            target_grid = ground.PARA.target_grid;
                        
            [dummy, pos_target] = min(abs(target_grid-top_pos));
            
            target_grid = target_grid(pos_target:end,1) - target_grid(pos_target,1);
            i=2;
            while target_grid(i,1) < sum(ground.STATVAR.layerThick,1) - 0.5 .* (target_grid(i,1)-target_grid(i-1,1))
                target_grid_new = target_grid(1:i,1);
                i=i+1;
                if i > size(target_grid,1)
                    target_grid = [target_grid; target_grid(end,1) + (target_grid(end,1) - target_grid(end-1,1))];
                end
            end
            target_grid = target_grid(1:i,1);
            target_grid(end,1) = sum(ground.STATVAR.layerThick,1);
            
            for i =1:size(variable_list,1)
                STATVAR_new.(variable_list{i,1}) = target_grid(1:end-1,1).*0;
            end
            
            %vector a original grid STATVAR.layerThick
            %vector_b new grid target_grid
            vector_a = [0; cumsum(ground.STATVAR.layerThick)];
            vector_b = target_grid;
            
            i_a = 1;
            %find starting point of a
            while i_a <= size(vector_a,1) && vector_a(i_a,1) < vector_b(1,1)
                i_a = i_a+1;
            end
            
            for i=1:size(vector_b,1)-1
                while i_a < size(vector_a,1) && vector_a(i_a+1,1) < vector_b(i+1,1)
                    fraction_in = (min(vector_b(i+1,1),vector_a(i_a+1,1)) - max(vector_b(i,1), vector_a(i_a,1))) ./ (vector_a(i_a+1,1) - vector_a(i_a,1));
                    
                    for j=1:size(variable_list, 1)
                        STATVAR_new.(variable_list{j,1})(i,1) = STATVAR_new.(variable_list{j,1})(i,1) + fraction_in .* ground.STATVAR.(variable_list{j,1})(i_a,1);
                    end
                    i_a=i_a+1;
                end
                i_a = min(i_a, size(vector_a, 1)-1);
                fraction_in = (min(vector_b(i+1,1),vector_a(i_a+1,1)) - max(vector_b(i,1), vector_a(i_a,1))) ./ (vector_a(i_a+1,1) - vector_a(i_a,1));
                for j=1:size(variable_list, 1)
                    STATVAR_new.(variable_list{j,1})(i,1) = STATVAR_new.(variable_list{j,1})(i,1) + fraction_in .* ground.STATVAR.(variable_list{j,1})(i_a,1);
                end
            end
            for i =1:size(variable_list,1)
                ground.STATVAR.(variable_list{i,1})  = STATVAR_new.(variable_list{i,1});
            end
            
            if ~strcmp(class(ground.NEXT), 'Bottom')  %set top_pos of next class to new depth
                ground.NEXT.STATVAR.top_depth_rel2groundSurface = top_pos + sum(ground.STATVAR.layerThick, 1);
            end
            
        end
        
        %regridding routine for snow
        function [snow, regridded_yesNo] = regrid_snow(snow, extensive_variables, intensive_variables, intensive_scaling_variable)
            
            %replaces:  snow = modify_grid(snow)
            regridded_yesNo = 0;
            
            %reduce
            if sum(double(snow.STATVAR.ice < 0.5.*snow.PARA.swe_per_cell.*snow.STATVAR.area)) > 0
                regridded_yesNo = 1;
                i=1;
                while i<size(snow.STATVAR.layerThick,1)
                    if snow.STATVAR.ice(i,1) < 0.5.* snow.PARA.swe_per_cell.*snow.STATVAR.area(i,1)
                        snow = merge_cells_intensive(snow, i, i+1, intensive_variables, intensive_scaling_variable);
                        snow = merge_cells_extensive(snow, i, i+1, extensive_variables);
                        %rest is done by diagostic step, get_T_water
                    end
                    i=i+1;
                end
                %last cell, i = size(snow.STATVAR.layerThick,1)
                if i > 1 && snow.STATVAR.ice(end,1) < 0.5 .* snow.PARA.swe_per_cell.*snow.STATVAR.area(end,1)
                    snow = merge_cells_intensive(snow, i-1, i, intensive_variables, intensive_scaling_variable);
                    snow = merge_cells_extensive(snow, i-1, i, extensive_variables);
                    %rest is done by diagostic step, get_T_water
                end
            end
            
            if snow.STATVAR.ice(1) > 1.5.*snow.PARA.swe_per_cell.*snow.STATVAR.area(1)  %expand, check only first cell
               
                regridded_yesNo = 1;
                split_fraction = snow.STATVAR.ice(1) ./ (snow.PARA.swe_per_cell.*snow.STATVAR.area(1)); %e.g. 1.6
                split_fraction = (split_fraction-1)./split_fraction;
                snow = split_cell_intensive(snow, 1, intensive_variables);
                snow = split_cell_extensive(snow, 1, split_fraction, extensive_variables);
                %rest is done by diagostic step, get_T_water
            end  
        end
        
        
        %--------service functions---------------
        %service function to split cells for extensive variables
        function ground = split_cell_extensive(ground, pos, split_fraction, variable_list)
            
            for i=1:size(variable_list, 1)
                ground.STATVAR.(variable_list{i,1}) = [ground.STATVAR.(variable_list{i,1})(1:pos-1,1); split_fraction .* ground.STATVAR.(variable_list{i,1})(pos,1); ...
                    (1-split_fraction) .* ground.STATVAR.(variable_list{i,1})(pos,1); ground.STATVAR.(variable_list{i,1})(pos+1:end,1)];
            end
        end
        
        %service function to split cells for intnsive variables
        function ground = split_cell_intensive(ground, pos, variable_list)
            
            for i=1:size(variable_list, 1)
                ground.STATVAR.(variable_list{i,1}) = [ground.STATVAR.(variable_list{i,1})(1:pos-1,1); ground.STATVAR.(variable_list{i,1})(pos,1); ...
                    ground.STATVAR.(variable_list{i,1})(pos,1); ground.STATVAR.(variable_list{i,1})(pos+1:end,1)];
            end
        end
        
        %service function to merge cells for extensive variables
        function ground = merge_cells_extensive(ground, pos1, pos2, variable_list)  %pos2 is deleted and merged with pos1; restriction: pos1 must be smaller than pos2
            
            for i=1:size(variable_list, 1)
                dummy = ground.STATVAR.(variable_list{i,1})(pos2,1);
                ground.STATVAR.(variable_list{i,1})(pos2,:) = [];
                ground.STATVAR.(variable_list{i,1}) = [ground.STATVAR.(variable_list{i,1})(1:pos1-1,1); ground.STATVAR.(variable_list{i,1})(pos1,1) + dummy; ...
                    ground.STATVAR.(variable_list{i,1})(pos1+1:end,1)];
            end
        end
        
        %service function to merge cells for intensive variables
        function ground = merge_cells_intensive(ground, pos1, pos2, variable_list, scaling_variable)  %pos2 is deleted and merged with pos1; restriction: pos1 must be smaller than pos2
            %weighted average, using scaling_variable as weighing function -> example: 
            %snow.STATVAR.d(i+1) = (snow.STATVAR.d(i+1).*snow.STATVAR.ice(i+1) + snow.STATVAR.d(i).*snow.STATVAR.ice(i))./(snow.STATVAR.ice(i+1) + snow.STATVAR.ice(i));
            for i=1:size(variable_list, 1)
                merged = (ground.STATVAR.(variable_list{i,1})(pos1,1) .* ground.STATVAR.(scaling_variable)(pos1,1) + ...
                    ground.STATVAR.(variable_list{i,1})(pos2,1) .* ground.STATVAR.(scaling_variable)(pos2,1)) ./ ...
                    (ground.STATVAR.(scaling_variable)(pos1,1) + ground.STATVAR.(scaling_variable)(pos2,1));
                ground.STATVAR.(variable_list{i,1})(pos2,:) = [];
                ground.STATVAR.(variable_list{i,1}) = [ground.STATVAR.(variable_list{i,1})(1:pos1-1,1); merged; ground.STATVAR.(variable_list{i,1})(pos1+1:end,1)];
            end
        end
        
        %same as above, but merging cells from two different objects
        function ground = merge_cells_extensive2(ground, pos, ground2, pos2, variable_list)  %cell at pos2 of ground2 is merged with pos of ground
            
            for i=1:size(variable_list, 1)
                ground.STATVAR.(variable_list{i,1})(pos,1) = ground.STATVAR.(variable_list{i,1})(pos,1) + ground2.STATVAR.(variable_list{i,1})(pos2,1);
            end
        end

        %same as above, but merging cells from two different objects        
        function ground = merge_cells_intensive2(ground, pos, ground2, pos2, variable_list, scaling_variable)  %cell at pos2 of ground2 is merged with pos of ground
            for i=1:size(variable_list, 1)
                ground.STATVAR.(variable_list{i,1})(pos,1) = (ground.STATVAR.(variable_list{i,1})(pos,1) .* ground.STATVAR.(scaling_variable)(pos,1) + ...
                    ground2.STATVAR.(variable_list{i,1})(pos2,1) .* ground2.STATVAR.(scaling_variable)(pos2,1)) ./ ...
                    (ground.STATVAR.(scaling_variable)(pos,1) + ground2.STATVAR.(scaling_variable)(pos2,1));
            end
        end
        
    end
end

    