%========================================================================
% CryoGrid class
% S. Westermann December 2021
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
        

    end
end