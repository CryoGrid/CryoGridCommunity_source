%========================================================================
% CryoGrid FORCING class FORCING_seb
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
%========================================================================

classdef FORCING_ubT < matlab.mixin.Copyable
    
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
          
            %temp=load(['forcing/' forcing.PARA.filename], 'FORCING');
            temp=load([forcing.PARA.forcing_path forcing.PARA.filename], 'FORCING');
            
            forcing.DATA.T_ub = temp.FORCING.data.T_ub;
            forcing.DATA.snow_depth = temp.FORCING.data.snow_depth;
            forcing.DATA.timeForcing = temp.FORCING.data.t_span;
            

            
            if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1)) >= 1e-10 %~=0
                disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
                forcing.STATUS=0;
                return
            else
                forcing.STATUS=1;
            end

            if isempty(forcing.PARA.start_time) || isnan(forcing.PARA.start_time(1,1)) %|| ~ischar(forcing.PARA.start_time)
                forcing.PARA.start_time = forcing.DATA.timeForcing(1,1);
            else
                %forcing.PARA.start_time = datenum(forcing.PARA.start_time, 'dd.mm.yyyy');
                forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            end
             
            if isempty(forcing.PARA.end_time) || isnan(forcing.PARA.end_time(1,1)) %~ischar(forcing.PARA.end_time)
                forcing.PARA.end_time = floor(forcing.DATA.timeForcing(end,1));
            else
                %forcing.PARA.end_time = datenum(forcing.PARA.end_time, 'dd.mm.yyyy');
                forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            end
            
            %initialize TEMP
            forcing.TEMP.T_ub=0;
            forcing.TEMP.snow_depth=0;

        end
        
        function forcing = interpolate_forcing(forcing, tile)
            t = tile.t;

            posit=floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;
            
            forcing.TEMP.snow_depth=forcing.DATA.snow_depth(posit,1)+(forcing.DATA.snow_depth(posit+1,1)-forcing.DATA.snow_depth(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.T_ub=forcing.DATA.T_ub(posit,1)+(forcing.DATA.T_ub(posit+1,1)-forcing.DATA.T_ub(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.T_ub = double(forcing.TEMP.snow_depth == 0 || forcing.TEMP.T_ub <0 ) .* forcing.TEMP.T_ub;
            
            forcing.TEMP.t = t;
        end
        


        
        
 
                
    end
end