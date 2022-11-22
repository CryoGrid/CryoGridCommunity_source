%========================================================================
% CryoGrid FORCING class FORCING_seb

% S. Westermann, November 2020
%========================================================================

classdef FORCING_read_out < matlab.mixin.Copyable
    
    properties
        forcing_index
        DATA            % forcing data time series
        TEMP            % forcing data interpolated to a timestep
        PARA            % parameters
        STATUS         
    end
    
    
    methods
        

            
        
        function forcing = provide_PARA(forcing)         

            %forcing.PARA.filename = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
        end
        
        function forcing = provide_CONST(forcing)
            
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end

        function forcing = finalize_init(forcing, tile)
            
            forcing.PARA.start_time = datenum(forcing.PARA.start_time, 'dd.mm.yyyy');
            
            forcing.PARA.end_time = datenum(forcing.PARA.end_time, 'dd.mm.yyyy');
        end


        function forcing = interpolate_forcing(forcing, tile)
        
        end

    end
end