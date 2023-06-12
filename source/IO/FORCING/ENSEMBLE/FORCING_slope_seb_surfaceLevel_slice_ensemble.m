%========================================================================
% CryoGrid FORCING class FORCING_seb_mat
%
% simple model forcing for GROUND classes computing the surface energy 
% balance (keyword "seb"). 
% 
% The data is obtained using the READ_FORCING_mat class. See this class for
% instructions about mat-file data format.
% 
% The mandatory forcing variables are:
%
% Tair:      Air temperature (in degree Celsius)
% Lin:       incoming long-wave radiation (in W/m2)
% Sin:       incoming short-wave radiation (in W/m2)
% rainfall:  Rainfall (in mm/day)
% snowfall:  Snowfall (in mm/day)
% q:         absolute humidity (in kg water vapor / kg air)
% p:         air pressure (OPTIONAL, in Pa)
% wind:       wind speed (in m/sec)
% 
% All forcing variables must be discretized identically, and one array of
% timestamps must be provided (t_span or timeForcing, in Matlab time / increment 1 
% corresponds to one day). 
%
% IMPORTANT POINT: the time series must be equally spaced in time, and this 
% must be really exact. When reading the timestamps from an existing data 
% set (e.g. an Excel file), rounding errors can result in small differences 
% in the forcing timestep, often less than a second off. In this case, it 
% is better to manually compile a new, equally spaced timestep in Matlab.
%
% Authors:
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
% T. Ingeman-Nielsen, October 2022
%
%========================================================================

classdef FORCING_slope_seb_surfaceLevel_slice_ensemble < FORCING_slope_seb_surfaceLevel_slice
    
    properties
        PERTURB_CLASS
    end
    
    methods
        
        function forcing = provide_PARA(forcing)         

            forcing = provide_PARA@FORCING_slope_seb_surfaceLevel_slice(forcing);
            
            forcing.PARA.perturb_forcing_class = [];   %perturbs forcing at runtime
            forcing.PARA.perturb_forcing_class_index = [];
        end
        
               
        function forcing = finalize_init(forcing, tile)
            forcing = finalize_init@FORCING_slope_seb_surfaceLevel_slice(forcing, tile);
            
            %optional post-processing with dedicated classes
            if ~isempty(forcing.PARA.perturb_forcing_class) && sum(isnan(forcing.PARA.perturb_forcing_class_index)) == 0
                forcing.PERTURB_CLASS = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(forcing.PARA.perturb_forcing_class){forcing.PARA.perturb_forcing_class_index,1});
                forcing.PERTURB_CLASS = finalize_init(forcing.PERTURB_CLASS, tile);
            end
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            forcing = interpolate_forcing@FORCING_slope_seb_surfaceLevel_slice(forcing, tile);
            forcing = perturb_forcing(forcing.PERTURB_CLASS, forcing);
        end


    end
end