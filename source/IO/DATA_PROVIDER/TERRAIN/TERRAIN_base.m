%========================================================================
% Base class for TERRAIN information, containing functions that can be 
% called from all inheriting, specific TERRAIN classes. 
% Thomas Ingeman-Nielsen, October 2022
%========================================================================

classdef TERRAIN_base < matlab.mixin.Copyable
    
    properties
        PARA
        CONST
        TEMP
        STATVAR
        DATA            % terrain data time series
    end
    
    
    methods
        
        function terrain = provide_PARA(terrain)
            terrain.PARA.altitude = [];     %
            terrain.PARA.slope_angle = [];     %
            terrain.PARA.aspect = [];     %
            terrain.PARA.skyview_factor = [];     %
            terrain.PARA.horizon_angles = [];     %
        end
        
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
        
        function terrain = initialize_terrain(terrain)
            % All angles in degrees and clockwise from N (maybe change)            
        end
        
    end
end