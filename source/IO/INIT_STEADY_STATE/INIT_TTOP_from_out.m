%========================================================================
% CryoGrid class
% S. Westermann December 2021
%========================================================================


classdef INIT_TTOP_from_out < matlab.mixin.Copyable
 

    properties
        PARA
        CONST
        STATVAR
	end
    
    
    methods

        
        function init = provide_PARA(init)         

            init.PARA.out_folder = []; %if empty, use the result folder
            init.PARA.out_file = [];

        end
        
        function init = provide_CONST(init)

        end
        
        function init = provide_STATVAR(init)

        end
		

		
		function init = finalize_init(init, tile)
             if isempty(init.PARA.out_folder) || sum(isnan(init.PARA.out_folder))>0
                 init.PARA.out_folder = [tile.PARA.result_path tile.PARA.run_name '/'];
                 init.PARA.out_file = [tile.PARA.run_name '_OUT_FDD_TDD.mat'];
             end
             temp=load([init.PARA.out_folder init.PARA.out_file], 'out');
             out = temp.out;
             tile.PARA.T_first_cell = out.STATVAR.TTOP;
             tile.PARA.start_depth_steady_state = out.STATVAR.TTOP_depth;
        end
        

    end
end