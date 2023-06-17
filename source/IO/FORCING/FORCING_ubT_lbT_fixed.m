%========================================================================
% CryoGrid FORCING class FORCING_ubT_lbT_fixed
%
% simple model forcing for GROUND classes depending only on fixed upper
% and lower boundary temperatures
%
% The mandatory userdefined parameters are:
% start_time
% end_time
% time_step
% T_ub   Upper boundary temperature.
% T_lb   Lower boundary temperature.
%
% T. Ingeman-Nielsen, December 2021
%========================================================================


classdef FORCING_ubT_lbT_fixed < FORCING_base
    
    properties
        forcing_index
        DATA            % forcing data time series
        TEMP            % forcing data interpolated to a timestep
        PARA            % parameters
        STATUS         
        CONST         
    end
    
    
    methods

        function forcing = provide_PARA(forcing)         
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.  

            forcing.PARA.start_time = [];     % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];       % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.time_step = [];      % time step of timeseries [days] (e.g. 0.5 days)
            forcing.PARA.T_ub = [];           % 
            forcing.PARA.T_lb = [];           % 
        end
        
        function forcing = provide_CONST(forcing)
            
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        function forcing = initialize_excel(forcing)
            
        end
 
        
        function forcing = finalize_init(forcing, tile)

            forcing = set_start_and_end_time(forcing); % assign start/end time

            forcing.STATUS = 1;

            % handle start time
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            
            % handle end time
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));

            % generate time sequence
            forcing.DATA.timeForcing = forcing.PARA.start_time:forcing.PARA.time_step:forcing.PARA.end_time;
            forcing.PARA.time_zero = datenum(datetime(year(datetime(forcing.PARA.start_time,'ConvertFrom','datenum')),1,1,0,0,0));

            %initialize TEMP
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Tair=0;
            forcing.TEMP.T_ub=0;
            forcing.TEMP.T_lb=0;
            
        end


        function forcing = interpolate_forcing(forcing, tile)
            % Interpolate forcing data to timestep tile.t
            t = tile.t;
            
            forcing.TEMP.T_ub = forcing.PARA.T_ub;
            forcing.TEMP.T_lb = forcing.PARA.T_lb;
            forcing.TEMP.Tair = forcing.TEMP.T_ub;
            forcing.TEMP.t = t;
        end

        function fig = plot(forcing)
			error('This function is not implemented/updated for this specific class')
        end

    end

end