% inhertits from SNOW_base_class and adds SEB as upper boundary
% designed to function as CHILD of a GROUND class that is compatible with
% SNOW classes; compatible with interaction classes IA_SNOW_GROUND and IA_SNOW_GROUND_fcSimple_salt

classdef SNOW_simple_seb_bucketW < SNOW_simple_seb

    methods

        %mandatory functions for each class
        
        function xls_out = write_excel(snow)
            xls_out = {'CLASS','index',NaN,NaN,NaN;'SNOW_simple_seb_bucketW',1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN;NaN,'value','default','unit',NaN;'density',300,350,'[kg/m3]','snow density';'albedo_max',0.850000000000000,0.850000000000000,'[-]','not active';'albedo_min',0.550000000000000,0.550000000000000,'[-]','not active';'albedo',0.800000000000000,0.800000000000000,'[-]','surface albedo';'epsilon',0.990000000000000,0.990000000000000,'[-]','surface emissivity';'z0',0.000100000000000000,0.000100000000000000,'[m]','roughness length';'field_capacity',0.0500000000000000,0.0500000000000000,'[-]','%fraction of porosity that can be filled with water before draining';'hydraulicConductivity',1,'    ','[m/sec]','    ';'dt_max',3600,3600,'[sec]','longest possible timestep';'dE_max',50000,50000,'[J/m3]','maximum change of energy per timestep';'swe_per_cell',0.0100000000000000,0.0100000000000000,'[m]','target SWE regulating grid cell size, 0.01m is ca. 3cm ';'CLASS_END',NaN,NaN,NaN,NaN};
        end

       function snow = provide_variables(snow)  %initializes the subvariables as empty arrays 
           snow = provide_variables@SNOW_simple_seb(snow); 
           snow = provide_PARA(snow);
           snow = provide_CONST(snow);
           snow = provide_STATVAR(snow);
        end
        
        function variable = initialize_from_file(snow, variable, section)
            variable = initialize_from_file@SNOW_simple_seb(snow, variable, section);
        end
        
        function snow = assign_global_variables(snow, forcing)
            snow = assign_global_variables@SNOW_simple_seb(snow, forcing);
        end
        
        function snow = initialize_zero_snow(snow, parentGround)
            snow = initialize_zero_snow@SNOW_simple_seb(snow, parentGround);
            %keyboard
            %snow.TEMP.d_water_ET = -snow.STATVAR.Qe ./ snow.CONST.L_s;
        end
        
        
        function snow = get_boundary_condition_u(snow, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            snow = get_boundary_condition_u@SNOW_simple_seb(snow, forcing);
            snow.TEMP.T_rainWater = forcing.TEMP.Tair;
            snow.TEMP.F_ub_water = snow.TEMP.rainfall;
            snow.TEMP.d_ice_sublim = -snow.STATVAR.Qe ./ snow.CONST.L_s; 
            snow.TEMP.d_E_sublim =  snow.TEMP.d_ice_sublim .* (snow.CONST.c_w .*snow.STATVAR.T(1) - snow.CONST.L_f);

       end
        
        function snow = get_boundary_condition_l(snow, forcing) 
            snow = get_boundary_condition_l@SNOW_simple_seb(snow, forcing);
            ground.TEMP.F_lb_water = 0; % zero flux if used as bottom class
        end
        
        
        function snow = get_derivatives_prognostic(snow)
            snow = get_derivatives_prognostic@SNOW_base_class(snow);
            snow = get_derivative_water(snow);      
            
        end
        
        function timestep = get_timestep(snow)  %could involve check for several state variables
             timestep = get_timestep@SNOW_simple_seb(snow);
        end
        
        function snow = advance_prognostic(snow, timestep) %real timestep derived as minimum of several classes
            snow = advance_prognostic@SNOW_base_class(snow, timestep);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* (snow.TEMP.snow_energy + snow.TEMP.d_E_sublim);  %rainfall energy already added
            snow.STATVAR.waterIce(1)  = snow.STATVAR.waterIce(1) + timestep .* (snow.TEMP.d_ice_sublim + snow.TEMP.snowfall) ; %rainfall done in advance_prognostic_water
            snow = advance_prognostic_water(snow, timestep);
            
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* (snow.TEMP.snowfall ./ (snow.PARA.density ./1000) + snow.TEMP.d_ice_sublim ./ (snow.STATVAR.ice(1) ./ snow.STATVAR.layerThick(1)));
            snow.STATVAR.target_density = min(1, snow.STATVAR.ice ./ snow.STATVAR.layerThick);
            snow.STATVAR.target_density(1) = min(1, (snow.STATVAR.ice(1) + timestep .* (snow.TEMP.snowfall + snow.TEMP.d_ice_sublim)) ./ snow.STATVAR.layerThick(1));
            
            if sum(snow.STATVAR.layerThick<=0)~=0 || sum(snow.STATVAR.waterIce<=0)~=0
                dff
            end
            
        end
        
        function snow = compute_diagnostic_first_cell(snow, forcing);
            snow = compute_diagnostic_first_cell@SNOW_simple_seb(snow, forcing);
        end
        
        function snow = compute_diagnostic(snow, forcing)
            snow = compute_diagnostic@SNOW_simple_seb(snow, forcing);
            
%             if sum(snow.STATVAR.layerThick==0)~=0 || sum(snow.STATVAR.waterIce<=0)~=0
%                  dff
%              end
            %snow = check_trigger@SNOW_simple_seb(snow);
        end
        
         function ground = troubleshoot(ground)
            ground = checkNaN(ground);
         end
            
        
        
        % non-mandatory
        function snow = get_derivative_water(snow)
            snow.TEMP.F_lb_water =0; %CHANGE LATER
            
            saturation = snow.STATVAR.water ./ max(1e-12, snow.STATVAR.layerThick - snow.STATVAR.ice);
            waterMobile = double(saturation > snow.PARA.field_capacity);
            snow.TEMP.d_water_out = waterMobile .* snow.PARA.hydraulicConductivity .* snow.STATVAR.water ./ snow.STATVAR.layerThick;
            snow.TEMP.d_water_in = snow.TEMP.d_water_out .*0;
            snow.TEMP.d_water_in(2:end,1) = snow.TEMP.d_water_out(1:end-1,1);
            snow.TEMP.d_water_in(1,1) = snow.TEMP.F_ub_water;
            if ~isempty(snow.IA_NEXT)
                get_boundary_condition_water_m(snow.IA_NEXT);
            else
                snow.TEMP.d_water_out(end,1) = snow.TEMP.F_lb_water; %CHANGE LATER
            end
        end
        
        function snow = advance_prognostic_water(snow, timestep)
            %snow.STATVAR.water  = snow.STAVAR.water + snow.TEMP.d_water_ET .* timestep; %add this function when SEB is present

            snow.TEMP.d_water_in = snow.TEMP.d_water_in .* timestep;
            snow.TEMP.d_water_out = snow.TEMP.d_water_out .* timestep;
            %limit outflow to field capacity
            snow.TEMP.d_water_out  = min(snow.TEMP.d_water_out, max(0, snow.STATVAR.water - snow.PARA.field_capacity .* (snow.STATVAR.layerThick -  snow.STATVAR.ice)));
            snow.TEMP.d_water_in(2:end,1) = snow.TEMP.d_water_out(1:end-1,1);
            %limit inflow so that unity is not exceeded
            snow.TEMP.d_water_in = max(0, min(snow.TEMP.d_water_in, snow.STATVAR.layerThick  - snow.STATVAR.waterIce));
            snow.TEMP.d_water_out(1:end-1,1) = snow.TEMP.d_water_in(2:end,1);
            
%            finalize_boundary_condition_water_m(snow.IA_NEXT);
            
            energy_out = snow.CONST.c_w .* snow.STATVAR.T .* snow.TEMP.d_water_out;
            energy_in = energy_out.*0;
            energy_in (2:end,1) = energy_out(1:end-1,1);
            energy_in(1,1) = snow.CONST.c_w .* snow.TEMP.T_rainWater .* snow.TEMP.d_water_in(1,1);
            
            snow.STATVAR.waterIce = snow.STATVAR.waterIce - snow.TEMP.d_water_out + snow.TEMP.d_water_in;
            snow.STATVAR.energy = snow.STATVAR.energy - energy_out + energy_in;

        end


    end
end