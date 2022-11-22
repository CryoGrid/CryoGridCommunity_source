%========================================================================
% TERRAIN class for cases where there is no local terrain (flat surface).
%
% Robin B. Zweigel, October 2022
%========================================================================

classdef TERRAIN_none < TERRAIN_base
    
    methods
        
        % function terrain = provide_PARA(terrain)
            % provide_PARA is inherited from base clase
        % end
        
        function terrain = provide_CONST(terrain)
            % empty; CONSTs are provided in the relevant terrain classes
        end
        
        function terrain = provide_STATVAR(terrain)
            % empty; STATVARs are provided in the relevant terrain classes
        end
        
        %%%%%% ----- finalize init functions ----- %%%%%%
        
        function terrain = initialize_TEMP(terrain)
            % The TEMP variables required for all classes
        end
        
        function terrain = finalize_init(terrain, tile)
            % All angles in degrees and clockwise from N (maybe change) 
            terrain.PARA.altitude = tile.PARA.altitude;     %
            terrain.PARA.slope_angle = 0;     %
            terrain.PARA.aspect = 0;     %
            terrain.PARA.skyview_factor = 1;     %
            terrain.PARA.horizon_angles = 0;     %
        end
        
    end
end