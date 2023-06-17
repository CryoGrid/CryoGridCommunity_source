%========================================================================
% CryoGrid FORCING class FORCING_ubT_mat
%
% simple model forcing for GROUND classes depending only on an upper
% boundary temperature (T_ub) forcing.
%
% The data is obtained using the READ_FORCING_mat class. See this class for
% instructions about mat-file data format.
%
% The mandatory forcing variables are:
%
% T_ub:        Upper boundary temperature (in degree Celsius)
%
% One array of timestamps must be provided (t_stamp, in Matlab time / 
% increment 1 corresponds to one day). 
%
% IMPORTANT POINT: the time series must be equally spaced in time, and this 
% must be really exact. When reading the timestamps from an existing data 
% set (e.g. an Excel file), rounding errors can result in small differences 
% in the forcing timestep, often less than a second off. In this case, it 
% is better to manually compile a new, equally spaced timestep in Matlab.
%
% T. Ingeman-Nielsen, October 2022
%========================================================================

classdef FORCING_ubT < FORCING_base & READ_FORCING_mat
        
    methods
        
        function forcing = provide_PARA(forcing)         

            forcing.PARA.filename = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain

        end

        function forcing = provide_CONST(forcing)
            
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        
        function forcing = finalize_init(forcing, tile)
          
            variables = {'T_ub'};
            [data, times] = read_mat([forcing.PARA.forcing_path forcing.PARA.filename], variables);
            
            for i=1:size(variables,1)
                if isfield(data, variables{i,1})
                    forcing.DATA.(variables{i,1}) = data.(variables{i,1});
                end
            end

            forcing = check_and_correct(forcing); % Remove known errors
            forcing = set_start_and_end_time(forcing); % assign start/end time
           
            %initialize TEMP
            forcing.TEMP.T_ub=0;
            forcing.TEMP.snow_depth=0;

        end
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_base(forcing, tile);
        end
        

        %-------------param file generation-----
        function forcing = param_file_info(forcing)
            forcing = provide_PARA(forcing);

            forcing.PARA.STATVAR = [];
            forcing.PARA.class_category = 'FORCING';
            
            forcing.PARA.comment.filename = {'filename of Matlab file containing forcing data'};
            
            forcing.PARA.default_value.forcing_path = {'forcing/'};
            forcing.PARA.comment.forcing_path = {'path where forcing data file is located'};
            
            forcing.PARA.comment.start_time = {'start time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.start_time.name =  'H_LIST';
            forcing.PARA.options.start_time.entries_x = {'year' 'month' 'day'};
            
            forcing.PARA.comment.end_time = {'end_time time of the simulations (must be within the range of data in forcing file) - year month day'};
            forcing.PARA.options.end_time.name =  'H_LIST'; % 
            forcing.PARA.options.end_time.entries_x = {'year' 'month' 'day'};
           
            forcing.PARA.default_value.heatFlux_lb = {0.05};
            forcing.PARA.comment.heatFlux_lb = {'heat flux at the lower boundary [W/m2] - positive values correspond to energy gain'};
            
        end

        
        
 
                
    end
end