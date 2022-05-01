%========================================================================
% CryoGrid INIT_STEADY_STATE class INIT_TTOP_from_forcing
% CryoGrid INIT_STEADY_STATE class used in conjcntion with the TILE_BUILDER
% class new_init_steady_state; calculates the temperature on the top of 
% permafrost (TTOP) from the model forcing (air temperatures) with the 
% TTOP approach. It is used in the accelerated spin-up procedure.

% S. Westermann, Jan 2022
%========================================================================


classdef INIT_TTOP_from_forcing < matlab.mixin.Copyable
 

    properties
        PARA
        CONST
        STATVAR
	end
    
    
    methods

        
        function init = provide_PARA(init)         

            init.PARA.nf = [];
            init.PARA.nt = [];
            init.PARA.rk = [];
            init.PARA.start_time = [];
            init.PARA.end_time = [];
        end
        
        function init = provide_CONST(init)

        end
        
        function init = provide_STATVAR(init)

        end
		

		
		function init = finalize_init(init, tile)
            tile.PARA.T_first_cell = [];
            
            %assign start and end time, user is responsible to select
            %meaningful dates
            if isempty(init.PARA.start_time) || isnan(init.PARA.start_time(1,1)) %|| ~ischar(forcing.PARA.start_time)
                init.PARA.start_time = tile.FORCING.PARA.start_time;
            else
                init.PARA.start_time = datenum(init.PARA.start_time(1,1), init.PARA.start_time(2,1), init.PARA.start_time(3,1));
            end
             
            if isempty(init.PARA.end_time) || isnan(init.PARA.end_time(1,1)) %~ischar(forcing.PARA.end_time)
                init.PARA.end_time = tile.FORCING.PARA.end_time;
            else
                init.PARA.end_time = datenum(init.PARA.end_time(1,1), init.PARA.end_time(2,1), init.PARA.end_time(3,1));
            end
            
            start_index = find( tile.FORCING.DATA.timeForcing(:,1) >= init.PARA.start_time, 1);
            end_index = find( tile.FORCING.DATA.timeForcing(:,1) >= init.PARA.end_time, 1);
            
            FDD = sum(tile.FORCING.DATA.Tair(start_index:end_index,1) .* double(tile.FORCING.DATA.Tair(start_index:end_index,1)<0),1);
            TDD = sum(tile.FORCING.DATA.Tair(start_index:end_index,1) .* double(tile.FORCING.DATA.Tair(start_index:end_index,1)>0),1);
            number_of_values = end_index - start_index + 1;
            
            
            is_PF = init.PARA.nt .* init.PARA.rk .* TDD + init.PARA.nf  .* FDD < 0;
            tile.PARA.T_first_cell = double(is_PF) .* (init.PARA.nt .* init.PARA.rk .* TDD + init.PARA.nf  .* FDD) ./ number_of_values + ...
                    double(~is_PF) .* (init.PARA.nt .* TDD + init.PARA.nf ./ init.PARA.rk  .* FDD) ./ number_of_values;

        end
        
        %-------------param file generation-----
        function init = param_file_info(init)
            init = provide_PARA(init);

            init.PARA.STATVAR = [];
            init.PARA.class_category = 'init_steady_state';

            init.PARA.default_value.nf = {0.5};
            init.PARA.comment.nf = {'winter n-factor'};
            
            init.PARA.default_value.nt = {1};  
            init.PARA.comment.nt = {'summer n-factor'};

            init.PARA.default_value.rk = {0.8};
            init.PARA.comment.rk = {'thawed divided by frozen thermal conductivity active layer'};         

            
            init.PARA.comment.start_time = {'start time of initialization period (must be within the range of data in forcing file) - year month day'};
            init.PARA.options.start_time.name =  'H_LIST';
            init.PARA.options.start_time.entries_x = {'year' 'month' 'day'};
            
            init.PARA.comment.end_time = {'end_time time of initialization period (must be within the range of data in forcing file) - year month day'};
            init.PARA.options.end_time.name =  'H_LIST'; % 
            init.PARA.options.end_time.entries_x = {'year' 'month' 'day'};

        end

    end
end