%========================================================================
% CryoGrid FORCING class FORCING_seb_ncFromTopoScale
%reads single nc-file produced by Python version of TopoScale
%========================================================================

classdef FORCING_seb_nc_PCCH < matlab.mixin.Copyable
    
    properties
        forcing_index
        DATA            % forcing data time series
        TEMP            % forcing data interpolated to a timestep
        PARA            % parameters
        STATUS
        CONST
    end
    
    
    methods
        
        %mandatory functions
        
        function forcing = provide_PARA(forcing)
            % INITIALIZE_PARA  Initializes PARA structure, setting the variables in PARA.
            forcing.PARA.filename = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
%             forcing.PARA.all_rain_T = [];
%             forcing.PARA.all_snow_T = [];
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
        end
        
        function forcing = provide_CONST(forcing)

        end
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
        
        function forcing = finalize_init(forcing, tile)
            
            
            %variables = {'t2m'; 'd2m'; 'u10'; 'v10'; 'ssrd'; 'strd'; 'tp'};
            variables = {'time'; 'Tair'; 'wind'; 'q'; 'Sin'; 'Lin'; 'p'; 'snowfall'; 'rainfall'};
            for i=1:size(variables,1)
                temp.(variables{i,1}) = double(squeeze(ncread([forcing.PARA.forcing_path forcing.PARA.filename], variables{i,1})));
            end
            temp.reference_time = ncinfo([forcing.PARA.forcing_path forcing.PARA.filename]);
            temp.reference_time =temp.reference_time.Variables(1).Attributes(3).Value(12:12+9);
%             temp.reference_time=temp.reference_time.Attributes(1).Value;
%             temp.reference_time =
%             datenum(temp.reference_time(end-18:end), 'yyyy-mm-dd
%             HH:MM:SS'); Still a problem in Simon's files, format seems to
%             change
            temp.reference_time = datenum(temp.reference_time, 'yyyy-mm-dd');
           
            temp.timeForcing = temp.reference_time + temp.time - temp.time(1) +0.25;
            
            temp.Tair = temp.Tair-273.15;
            temp.snowfall = temp.snowfall * 8;
            temp.rainfall = temp.rainfall * 8;
            temp.Lin = temp.Lin./3600./3;
            temp.Sin = temp.Sin./3600./3;
            temp.time = [];
            
            forcing.DATA = temp;
            
            if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1))>1e-9
                disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
                forcing.STATUS=0;
                return
            else
                forcing.STATUS=1;
            end
            
            forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
            forcing.DATA.Lin(find(forcing.DATA.Lin==0)) = 5.67e-8 .* (forcing.DATA.Tair(find(forcing.DATA.Lin==0))+273.15).^4;
            forcing.DATA.Sin(forcing.DATA.Sin<0)=0;
            forcing.DATA.rainfall=forcing.DATA.rainfall.*forcing.PARA.rain_fraction;
            forcing.DATA.snowfall=forcing.DATA.snowfall.*forcing.PARA.snow_fraction;
            
            if isempty(forcing.PARA.start_time) || isnan(forcing.PARA.start_time(1,1)) % || ~ischar(forcing.PARA.start_time)
                forcing.PARA.start_time = forcing.DATA.timeForcing(1,1);
            else
                forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            end
            if isempty(forcing.PARA.end_time) || isnan(forcing.PARA.end_time(1,1)) %|| ~ischar(forcing.PARA.end_time)
                forcing.PARA.end_time = floor(forcing.DATA.timeForcing(end,1));
            else
                forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1),forcing.PARA.end_time(3,1));
            end
            

            
            %initialize TEMP
            forcing.TEMP.snowfall=0;
            forcing.TEMP.rainfall=0;
            forcing.TEMP.Lin=0;
            forcing.TEMP.Sin=0;
            forcing.TEMP.Tair=0;
            forcing.TEMP.wind=0;
            forcing.TEMP.q=0;
            forcing.TEMP.p=0;
            
        end
        
        
        function forcing = interpolate_forcing(forcing, tile)
            t = tile.t;
            
            posit=floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;
            
            forcing.TEMP.snowfall=forcing.DATA.snowfall(posit,1)+(forcing.DATA.snowfall(posit+1,1)-forcing.DATA.snowfall(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.rainfall=forcing.DATA.rainfall(posit,1)+(forcing.DATA.rainfall(posit+1,1)-forcing.DATA.rainfall(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.Lin=forcing.DATA.Lin(posit,1)+(forcing.DATA.Lin(posit+1,1)-forcing.DATA.Lin(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.Sin=forcing.DATA.Sin(posit,1)+(forcing.DATA.Sin(posit+1,1)-forcing.DATA.Sin(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.Tair=forcing.DATA.Tair(posit,1)+(forcing.DATA.Tair(posit+1,1)-forcing.DATA.Tair(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.wind=forcing.DATA.wind(posit,1)+(forcing.DATA.wind(posit+1,1)-forcing.DATA.wind(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.q=forcing.DATA.q(posit,1)+(forcing.DATA.q(posit+1,1)-forcing.DATA.q(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.p=forcing.DATA.p(posit,1)+(forcing.DATA.p(posit+1,1)-forcing.DATA.p(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
            
            forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
            forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
            forcing.TEMP.t = t;
            
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
            
            forcing.PARA.default_value.rain_fraction = {1};  
            forcing.PARA.comment.rain_fraction = {'rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)'};
            
            forcing.PARA.default_value.snow_fraction = {1};  
            forcing.PARA.comment.snow_fraction = {'snowfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)'};

%             forcing.PARA.default_value.all_rain_T = {0.5};
%             forcing.PARA.comment.all_rain_T = {'temperature above which all precip is rain'};
%                         
%             forcing.PARA.default_value.all_snow_T = {-0.5};
%             forcing.PARA.comment.all_snow_T = {'temperature below which all precip is snow'};
%             
%             forcing.PARA.default_value.slope_angle = {0}; 
%             forcing.PARA.comment.slope_angle = {'slope angle in degree'};
%             
%             forcing.PARA.default_value.aspect = {90}; 
%             forcing.PARA.comment.aspect = {'aspect of the slope in degrees'};
%             
%             forcing.PARA.default_value.albedo_surrounding_terrain = {0.2};
%             forcing.PARA.comment.albedo_surrounding_terrain = {'albedo of field of view from where solar radiation is reflected'};
%             
%             forcing.PARA.default_value.sky_view_factor ={1}; 
%             forcing.PARA.comment.sky_view_factor = {'sky view factor (0.5 for vertical rock walls)'};
            
            forcing.PARA.default_value.heatFlux_lb = {0.05};
            forcing.PARA.comment.heatFlux_lb = {'heat flux at the lower boundary [W/m2] - positive values correspond to energy gain'};
            
            forcing.PARA.default_value.airT_height = {2};  
            forcing.PARA.comment.airT_height = {'height above ground surface where air temperature from forcing data is applied'};

        end
        
    end
end