% base class build a model tile

classdef RUN_3D_Marco  < RUN_3D_STANDARD
    

    
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)


            run_info = provide_PARA@RUN_3D_STANDARD(run_info);
            run_info.PARA.altitudes =[];
            
        end
%         
%         function run_info = provide_CONST(run_info)
% 
%         end
%         
%         function run_info = provide_STATVAR(run_info)
% 
%         end
%         
%         function run_info = initialize_excel(run_info)
%             
%         end
%         
%         function run_info = finalize_init(run_info)
%  
%         end
%         
%         
        
%         function [run_info, tile] = run(run_info)
%             %this could first open spmd and assign run_number depending on
%             %worker, then do another round of pprovider
%             %it could also do a loop over different tile representing
%             %different sections of the run, e.g. initial inial init, spin-up, actual run 
%             parpool(run_info.PARA.number_of_tiles)
%             spmd
%                 run_info.PARA.worker_number = labindex;
%                 
%                 %update the worker-specific name of the parameter file and the run
%                 run_info.PPROVIDER = update_parameter_file(run_info.PPROVIDER, run_info.PARA.param_file_number(run_info.PARA.worker_number,1));
%                 run_info.PPROVIDER = update_run_name(run_info.PPROVIDER, run_info.PARA.worker_number);
%                 %read worker-specific parameter file
%                 run_info.PPROVIDER = read_parameters(run_info.PPROVIDER);
%                 
%                 run_info = customize(run_info);
%                 
%                 %tile = run_info.PPROVIDER.FUNCTIONAL_CLASSES.TILE{run_info.PARA.tile_number,1};
%                 tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
% 
%                 tile.RUN_INFO = run_info;
%                 run_info.TILE = tile;
% 
%                 tile = finalize_init(tile);
%                 
%                 tile = run(tile);  %time integration
%             end
%         end
 
        
        function run_info = customize(run_info)
            run_info.PPROVIDER.CLASSES.FORCING_seb_netcdf{1,1}.PARA.altitude = run_info.altitude(run_info.PARA.worker_number);

        end
        
        
        
        
        
    end
end



