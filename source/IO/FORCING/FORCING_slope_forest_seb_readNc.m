%========================================================================
% CryoGrid FORCING class FORCING_slope_seb
% simple model forcing for GROUND classes computing the surface energy balance 
% (keyword “seb”). The data must be stored in a Matlab “.mat” file which contains 
% a struct FORCING with field “data”, which contain the time series of the actual 
% forcing data, e.g. FORCING.data.Tair contains the time series of air temperatures. 
% Have a look at the existing forcing files in the folder “forcing” and prepare 
% new forcing files in the same way. The mandatory forcing variables are air temperature 
% (Tair, in degree Celsius), incoming long-wave radiation (Lin, in W/m2), 
% incoming short-.wave radiation (Sin, in W/m2), absolute humidity (q, in 
% kg water vapor / kg air), wind speed (wind, in m/sec), rainfall (rainfall, in mm/day), 
% snowfall (snowfall, in mm/day) and timestamp (t_span, 
% in Matlab time / increment 1 corresponds to one day). 
% IMPORTANT POINT: the time series must be equally spaced in time, and this must be 
% really exact. When reading the timestamps from an existing data set (e.g. an Excel file),
% rounding errors can result in small differences in the forcing timestep, often less 
% than a second off. In this case, it is better to manually compile a new, equally spaced 
% timestep in Matlab.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
% Edited (changes for slopes) were made by J. Schmidt, December 2020
% The forcing data need in addition for the slope:
%       S_TOA: Short-wave radiation at the top of the atmosphere
%       albedo_foot: albedo at the foot of the slope (can vary in time for
%           example with a snow layer
%   Not-mandatory (just applicable if there is a water body at the foot of
%   the slope):
%       seaT: seawater temperature
%       seaIce: 0 = time steps without sea ice; 1 = time steps with sea ice
%========================================================================

classdef FORCING_slope_forest_seb_readNc < FORCING_slope_seb_readNc %matlab.mixin.Copyable
    
    properties
      
    end
    
    
    methods
        
        %mandatory functions
        
        function forcing = provide_PARA(forcing)         

            forcing = provide_PARA@FORCING_slope_seb_readNc(forcing);
            
            % e.g. remote sensing data. Call this field "fc" (between 0 and 1).
% Also requires canopy transmissivity "tauc", canopy height "hc", and a 
% canopy extinction coefficient "muc". Example values are provided below.
% See Link & Marks (1999), Garren & Marks (2005), and Bair et al. (2016)
% for more on these parametrizations.
            
            forcing.PARA.canopy_transmissivity = [];
            forcing.PARA.canopy_height = []; % meters.
            forcing.PARA.canopy_extinction_coefficient = [];
            forcing.PARA.emissivity_canopy = []; % Canopy emissivity.
            forcing.PARA.fractional_canopy_covered_area = [];  % fractional canopy-covered area, should be read in from external data.

        end
        
        function forcing = provide_CONST(forcing)
            forcing.CONST.Tmfw = [];
            forcing.CONST.sigma = [];
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end

        
        function forcing = finalize_init(forcing, tile)
            

             forcing = finalize_init@FORCING_slope_seb_readNc(forcing, tile);
             
             % This outlines potential canopy corrections that could be done for the
             % radiation terms. We may want to include this in the downscaling routine
             % itself, obtaining one set of radiation fields for canopy-covered areas
             % another for open areas.
             
             % Requires a (static or dynamic) fractional canopy-covered area field from
             % e.g. remote sensing data. Call this field "fc" (between 0 and 1).
             % Also requires canopy transmissivity "tauc", canopy height "hc", and a
             % canopy extinction coefficient "muc". Example values are provided below.
             % See Link & Marks (1999), Garen & Marks (2005), and Bair et al. (2016)
             % for more on these parametrizations.
             
             % Canopy parameters (may have to read these in from external data). These
             % are applicable to Mammoth lakes.
             
             forcing.DATA.Lin = forcing.PARA.fractional_canopy_covered_area.*(forcing.PARA.canopy_transmissivity .* forcing.DATA.Lin + ...
                 (1 - forcing.PARA.canopy_transmissivity) .* forcing.PARA.emissivity_canopy .* forcing.CONST.sigma .* (forcing.DATA.Tair+forcing.CONST.Tmfw).^4) ...
                 + (1 - forcing.PARA.fractional_canopy_covered_area) .* forcing.DATA.Lin;
             
             sun_angle = forcing.DATA.sunElevation .* double(forcing.DATA.sunElevation > 0);
             forcing.DATA.Sin_dir = forcing.PARA.fractional_canopy_covered_area .* (forcing.DATA.Sin_dir .* ...
                 exp(-forcing.PARA.canopy_extinction_coefficient .* forcing.PARA.canopy_height ./ max(1e-12, sind(sun_angle)))) ...
                 + (1-forcing.PARA.fractional_canopy_covered_area) .* forcing.DATA.Sin_dir;

             forcing.DATA.Sin_dif = forcing.PARA.fractional_canopy_covered_area.*(forcing.DATA.Sin_dif.*forcing.PARA.canopy_transmissivity) + ...
                 (1 - forcing.PARA.fractional_canopy_covered_area).*forcing.DATA.Sin_dif;
             forcing.DATA.Sin = forcing.DATA.Sin_dif + forcing.DATA.Sin_dir;
             
        end
            

        function forcing = interpolate_forcing(forcing, tile)
            
            forcing = interpolate_forcing@FORCING_slope_seb_readNc(forcing, tile);
            
        end
        


                
        end
end