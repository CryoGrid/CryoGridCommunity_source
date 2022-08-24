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

classdef FORCING_seb < matlab.mixin.Copyable
    
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

            forcing.PARA.filename = [];   %filename of Matlab file containing forcing data
            forcing.PARA.forcing_path = [];
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            forcing.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
%             forcing.PARA.latitude = [];  % latitude
%             forcing.PARA.longitude = []; % longitude
%             forcing.PARA.altitude = [];  % elevation above sea level [m]
%             forcing.PARA.domain_depth = []; % total depth of the model domain [m]
            forcing.PARA.heatFlux_lb = [];  % heat flux at the lower boundary [W/m2] - positive values correspond to energy gain
            forcing.PARA.airT_height = [];  % height above ground at which air temperature (and wind speed!) from the forcing data are applied.
%             forcing.PARA.area = [];         % area of the grid cell [m2] 
            forcing.PARA.slope_angle = [];
        end
        
        
        
        function forcing = provide_CONST(forcing)
            
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end
        
%         function forcing = initialize_excel(forcing)
%             
%         end
        
       
        
        
        function forcing = finalize_init(forcing, tile)
          
            %temp=load(['forcing/' forcing.PARA.filename], 'FORCING');
            temp=load([forcing.PARA.forcing_path forcing.PARA.filename], 'FORCING');
            
            forcing.DATA.rainfall=temp.FORCING.data.rainfall.*forcing.PARA.rain_fraction;
            forcing.DATA.snowfall=temp.FORCING.data.snowfall.*forcing.PARA.snow_fraction;
            forcing.DATA.Tair = temp.FORCING.data.Tair;
            forcing.DATA.Lin = temp.FORCING.data.Lin;
            forcing.DATA.Sin = temp.FORCING.data.Sin;
            forcing.DATA.q = temp.FORCING.data.q;
            forcing.DATA.wind = temp.FORCING.data.wind;
            forcing.DATA.timeForcing = temp.FORCING.data.t_span;
            
            % Update spacial data if included in forcing file
%             if isfield(temp.FORCING.data,'z')
%                 forcing.PARA.altitude = round(temp.FORCING.data.z);
%             end
%             if isfield(temp.FORCING.data,'lon') &&isfield(temp.FORCING.data,'lat')
%                 forcing.PARA.longitude = temp.FORCING.data.lon;
%                 forcing.PARA.latitude  = temp.FORCING.data.lat;
%             end
            
            if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1)) >= 1e-10 %~=0
                disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
                forcing.STATUS=0;
                return
            else
                forcing.STATUS=1;
            end
            
            %here, consistency checks, RH->q calculation, set threhsolds for wind, etc. could be placed
            
            forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
            forcing.DATA.Lin(find(forcing.DATA.Lin==0)) = 5.67e-8 .* (forcing.DATA.Tair(find(forcing.DATA.Lin==0))+273.15).^4;
            
            %set pressure to mean pressure at corresponding altitude (international
            %altitude formula) if now provided by the forcing time series
            if ~isfield(temp.FORCING.data, 'p')
                forcing.DATA.p=forcing.DATA.Tair.*0 + 1013.25.*100.*(1-0.0065./288.15.*tile.PARA.altitude).^5.255;
            else
                forcing.DATA.p = temp.FORCING.data.p;
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
        

        function xls_out = write_excel(forcing)
			% XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
			
            xls_out = {'FORCING','index',NaN,NaN;'FORCING_seb',1,NaN,NaN;NaN,NaN,NaN,NaN;'filename',NaN,NaN,NaN;'start_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the first timestamp of the forcing data set will be used';'end_time',NaN,NaN,'provide in format dd.mm.yyyy; if left empty, the last timestamp of the forcing data set will be used';'rain_fraction',1,'[-]','rainfall in forcing file multiplied by this number';'snow_fraction',1,'[-]','snowfall in forcing file multiplied by this number';'latitude',NaN,'[degree]','geographical coordinates';'longitude',NaN,'[degree]',NaN;'altitude',NaN,'[m]','a.s.l.';'domain_depth',100,'[m]','should match a GRID point, model domain extends to this depth';'heatFlux_lb',0.0500000000000000,'[W/m2]','geothermal heat flux';'airT_height',2,'[m]','height of air temperature';'FORCING_END',NaN,NaN,NaN};
        end
        

        
                %constructor
%         function forcing = FORCING_seb(varargin)               % Temporary definition, to allow old code to run
%         %function forcing = FORCING_seb(index, pprovider, fprovider)      % Definition to be used when old code is no longer supported
%             % CONSTRUCTOR for FORCING_seb
%             %   Reads in forcing data from the specified file.
%             %
%             %   ARGUMENTS:
%             %   index:      user defined class index
%             %   pprovider:  instance of PARAMETER_PROVIDER class
% 			%   fprovider:  instance of FORCING_PROVIDER class
%             
%             % The following is only needed to allow legacy code to run
%             % May be removed when deprecated functions are removed
%             if nargin==3
%                 index = varargin{1};
%                 pprovider = varargin{2};
% 				fprovider = varargin{3};
%             else
%                 st = dbstack;
%                 warning(['DEPRECATION WARNING: Instantiating ' st.name '() with no arguments is deprecated.' newline,...
%                          'You should update your code to take advantage of new IO interface.']);
%                 return
%             end
%             % End allow legacy code
%             
%             forcing.forcing_index = index;
%             forcing = forcing.initialize();
%             forcing = forcing.populate_PARA(pprovider);
%             forcing = forcing.populate_DATA(fprovider);
%             forcing = forcing.finalize_setup();
%         end

%         
%         function forcing = initialize(forcing)
%             % INITIALIZE  Initializes all properties needed by the class.
% 
%             forcing.STATUS = 0;
%             forcing = forcing.initialize_PARA();
%             forcing = forcing.initialize_TEMP();
%             forcing = forcing.initialize_DATA();
%         end
%      
        
%         function forcing = populate_PARA(forcing, pprovider)
%             % POPULATE_PARA  Updates the PARA structure with values from pprovider.
%             %
%             %   ARGUMENTS:
%             %   pprovider:  instance of PARAMETER_PROVIDER class
%             
%             forcing.PARA = pprovider.populate_struct(forcing.PARA, 'FORCING', mfilename('class'), forcing.forcing_index);
%             if ~isfield(forcing.PARA, 'filepath')
%                 % if path to forcing data is not specified in PARA
%                 % assume it is in a forcing sub-directory to the main
%                 % script
%                 forcing.PARA.filepath = './forcing';
%             end
%         end
%         
% 
%         function forcing = populate_DATA(forcing, fprovider)
%             % POPULATE_DATA  Updates the DATA structure with values from fprovider.
%             %
%             %   ARGUMENTS:
% 			%	fprovider:	instance of FORCING class
% 
%             forcing.DATA = fprovider.populate_struct(forcing.DATA);
%             forcing.DATA.timeForcing = forcing.DATA.t_span;
%         end
            
        
%         function forcing = finalize_init(forcing)
%             % FINALIZE_SETUP  Performs all additional property
%             %   initializations and modifications. Checks for some (but not
%             %   all) data validity.
%             
%             % apply rain and snow fractions
%             forcing.DATA.rainfall = forcing.DATA.rainfall.*forcing.PARA.rain_fraction;
%             forcing.DATA.snowfall = forcing.DATA.snowfall.*forcing.PARA.snow_fraction;        
%             
%             % Update spatial data if included in forcing DATA
%             if ~isempty(forcing.DATA.z)
%                 forcing.PARA.altitude = round(forcing.DATA.z);
%             end
%             if ~isempty(forcing.DATA.lon) && ~isempty(forcing.DATA.lat)
%                  forcing.PARA.longitude = forcing.DATA.lon;
%                  forcing.PARA.latitude  = forcing.DATA.lat;
%             end
% 
%             if std(forcing.DATA.timeForcing(2:end,1)-forcing.DATA.timeForcing(1:end-1,1))>=1e-10
%                 error('timestamp of forcing data is not in regular intervals -> check, fix and restart')
%                 %forcing.STATUS=0;
%                 %return
%             %else
%             %    forcing.STATUS=1;
%             end
%             
%             forcing.STATUS = 1;
% 
%             % here, consistency checks, RH->q calculation, set threhsolds for wind, etc. could be placed
% 
%             % forcing.DATA.wind(forcing.DATA.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence
%             forcing.DATA.Lin(find(forcing.DATA.Lin==0)) = 5.67e-8 .* (forcing.DATA.Tair(find(forcing.DATA.Lin==0))+273.15).^4;
% 
%             %set pressure to mean pressure at corresponding altitude (international
%             %altitude formula) if now provided by the forcing time series
%             if isempty(forcing.DATA.p)
%                 forcing.DATA.p = forcing.DATA.Tair.*0 + 1013.25.*100.*(1-0.0065./288.15.*forcing.PARA.altitude).^5.255;
%             end
% 
%             if isempty(forcing.PARA.start_time) || ~ischar(forcing.PARA.start_time)
%                 forcing.PARA.start_time = forcing.DATA.timeForcing(1,1);
%             else
%                 forcing.PARA.start_time = datenum(forcing.PARA.start_time, 'dd.mm.yyyy');
%             end
%             if isempty(forcing.PARA.end_time) || ~ischar(forcing.PARA.end_time)
%                 forcing.PARA.end_time = floor(forcing.DATA.timeForcing(end,1));
%             else
%                 forcing.PARA.end_time = datenum(forcing.PARA.end_time, 'dd.mm.yyyy');
%             end
%         end
        
        
 
                
    end
end