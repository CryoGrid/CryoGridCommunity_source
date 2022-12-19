%========================================================================
% CryoGrid FORCING class FORCING_do_nothing
% simple model forcing for GROUND classes computing the surface energy 
% balance  (keyword "seb"). 
% All parameters are hardcoded to 0, except pressure which is hardcoded
% 1e5.
%
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================

classdef FORCING_do_nothing < FORCING_base
    
    methods
        
        function forcing = provide_PARA(forcing)         
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
        end
        

        function forcing = provide_CONST(forcing)
            
        end
        
        
        function forcing = provide_STATVAR(forcing)

        end
        
        
        function forcing = finalize_init(forcing, tile)
            
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));

        end
        

        function forcing = interpolate_forcing(forcing, tile)
            t = tile.t;
            
            forcing.TEMP.snowfall=0;           
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Lin=0;             
            forcing.TEMP.Sin=0 ;
            forcing.TEMP.Tair=0;
            forcing.TEMP.wind=0;
            forcing.TEMP.q=0;
            forcing.TEMP.p=1e5;
            forcing.TEMP.t = t;
        end

        
        %-------------param file generation-----
        function forcing = param_file_info(forcing)
            forcing = provide_PARA(forcing);

            forcing.PARA.STATVAR = [];
            forcing.PARA.class_category = 'FORCING';
            forcing.PARA.default_value = [];
            
            forcing.PARA.comment.start_time = {'start time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.start_time.name =  'H_LIST';
            forcing.PARA.options.start_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.comment.end_time = {'end_time time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.end_time.name =  'H_LIST'; % 
            forcing.PARA.options.end_time.entries_x = {'year' 'month' 'day'};
            
        end
                
    end
end