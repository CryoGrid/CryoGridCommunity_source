% inhertits from SNOW_base_class and adds SEB as upper boundary
% designed to function as CHILD of a GROUND class that is compatible with
% SNOW classes; compatible with interaction classes IA_SNOW_GROUND and IA_SNOW_GROUND_fcSimple_salt

classdef SNOW_simple_seb < SNOW_base_class

    methods

        %mandatory functions for each class
        function xls_out = write_excel(snow)
            xls_out = {'CLASS','index',NaN,NaN,NaN;'SNOW_simple_seb',1,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN;NaN,'value','default','unit',NaN;'density',300,350,'[kg/m3]','snow density';'albedo_max',0.850000000000000,0.850000000000000,'[-]','not active';'albedo_min',0.550000000000000,0.550000000000000,'[-]','not active';'albedo',0.800000000000000,0.800000000000000,'[-]','surface albedo';'epsilon',0.990000000000000,0.990000000000000,'[-]','surface emissivity';'z0',0.000100000000000000,0.000100000000000000,'[m]','roughness length';'field_capacity',0.0500000000000000,0.0500000000000000,'[-]','%fraction of porosity that can be filled with water before draining';'hydraulicConductivity',1,'    ','[m/sec]','    ';'dt_max',3600,3600,'[sec]','longest possible timestep';'dE_max',50000,50000,'[J/m3]','maximum change of energy per timestep';'swe_per_cell',0.0100000000000000,0.0100000000000000,'[m]','target SWE regulating grid cell size, 0.01m is ca. 3cm ';'CLASS_END',NaN,NaN,NaN,NaN};
        end


       function snow = provide_variables(snow)  %initializes the subvariables as empty arrays 
           snow = provide_variables@SNOW_base_class(snow); 
           snow = provide_PARA(snow);
           snow = provide_CONST(snow);
           snow = provide_STATVAR(snow);
        end
        
        function variable = initialize_from_file(snow, variable, section)
            variable = initialize_from_file@SNOW_base_class(snow, variable, section);
        end
        
        function snow = assign_global_variables(snow, forcing)
            snow.PARA.airT_height = forcing.PARA.airT_height;
        end
        
        function snow = initialize_zero_snow(snow, parentGround)
            snow = initialize_zero_snow@SNOW_base_class(snow, parentGround);
        end
        
        
        function snow = get_boundary_condition_u(snow, forcing) %functions specific for individual class, allow changing from Dirichlet to SEB
            snow = surface_energy_balance(snow, forcing);
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000./(24.*3600); % Simone!! ./(24.*3600); %snowfall is in mm/day
            snow.TEMP.rainfall = forcing.TEMP.rainfall ./1000./(24.*3600); % Simone!! ./(24.*3600);
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);
            snow.TEMP.rain_energy = snow.TEMP.rainfall .* max(0, forcing.TEMP.Tair) .* snow.CONST.c_w;
        end
        
        function snow = get_boundary_condition_u_CHILD(snow, forcing)
             snow = get_boundary_condition_u(snow, forcing); %same function as  for normal snow class
        end
        
        function snow = get_boundary_condition_u_create_CHILD(snow, forcing) 
            snow.TEMP.snowfall = forcing.TEMP.snowfall ./1000./(24.*3600); % Simone!! ./(24.*3600); %snowfall is in mm/day
            snow.TEMP.snow_energy = snow.TEMP.snowfall .* (min(0, forcing.TEMP.Tair) .* snow.CONST.c_i - snow.CONST.L_f);
        end
        
        function snow = get_boundary_condition_l(snow, forcing) 
            snow = get_boundary_condition_l@SNOW_base_class(snow, forcing);
        end
        
        function snow = get_derivatives_prognostic(snow)
            snow = get_derivatives_prognostic@SNOW_base_class(snow);
            snow.TEMP.d_energy(1) = snow.TEMP.d_energy(1) + snow.TEMP.rain_energy;  %add this here, since it can melt snow and does not change layerThick - must be taken into account for timestep calculation
        end
        
        function snow = get_derivatives_prognostic_CHILD(snow)  
            snow.TEMP.d_energy = snow.TEMP.F_ub + snow.TEMP.snow_energy + snow.TEMP.rain_energy;
        end
        
        function timestep = get_timestep(snow)  %could involve check for several state variables
             timestep1 = get_timestep@SNOW_base_class(snow);
             timestep2 = min((-snow.STATVAR.energy ./ snow.TEMP.d_energy) .*double(snow.TEMP.d_energy>0) + double(snow.TEMP.d_energy<=0).*1e5); %when snow is melting, do not melt more than there is in a grid cell
             timestep = min(timestep1, timestep2);
        end
        
        function timestep = get_timestep_CHILD(snow)  %will be ignored if it has a negative value
            %timestep = get_timestep(snow);
            timestep = -snow.STATVAR.energy ./ snow.TEMP.d_energy;
        end
        
        
        function snow = advance_prognostic(snow, timestep) %real timestep derived as minimum of several classes
            snow = advance_prognostic@SNOW_base_class(snow, timestep);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + timestep .* snow.TEMP.snow_energy;  %rainfall energy already added
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + timestep .* (snow.TEMP.snowfall + snow.TEMP.rainfall);
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./1000);
            snow.STATVAR.target_density = snow.STATVAR.ice ./ snow.STATVAR.layerThick;
            snow.STATVAR.target_density(1) = (snow.STATVAR.ice(1) + timestep .* snow.TEMP.snowfall) ./ snow.STATVAR.layerThick(1);
        end
        
        function snow = advance_prognostic_create_CHILD(snow, timestep)
            snow.STATVAR.energy = timestep .* snow.TEMP.snow_energy;
            snow.STATVAR.waterIce = timestep .* snow.TEMP.snowfall;
            snow.STATVAR.layerThick = timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./1000);
            snow.STATVAR.water = 0;
            snow.STATVAR.ice = snow.STATVAR.waterIce;
            snow.STATVAR.target_density = snow.STATVAR.ice ./ snow.STATVAR.layerThick;
        end
        
        function snow = advance_prognostic_CHILD(snow, timestep)
            
            snow.STATVAR.energy = snow.STATVAR.energy + timestep .* snow.TEMP.d_energy;
            snow.STATVAR.waterIce = snow.STATVAR.waterIce + timestep .* (snow.TEMP.snowfall + snow.TEMP.rainfall);
            snow.STATVAR.layerThick = snow.STATVAR.layerThick + timestep .* snow.TEMP.snowfall ./ (snow.PARA.density ./1000);
            snow.STATVAR.target_density = min(1,(snow.STATVAR.ice + timestep .* snow.TEMP.snowfall) ./ snow.STATVAR.layerThick);
      
        end
        
        function snow = compute_diagnostic_first_cell(snow, forcing)
            snow = L_star(snow, forcing);
        end
        
        function snow = compute_diagnostic(snow, forcing)
            snow = compute_diagnostic@SNOW_base_class(snow, forcing);
            snow = check_trigger(snow); 
        end
        
                
        function snow = compute_diagnostic_CHILD(snow, forcing)
            
            snow = get_T_water(snow);
            snow = conductivity(snow);
            
            snow.STATVAR.upperPos = snow.STATVAR.lowerPos + snow.STATVAR.layerThick;

            if snow.STATVAR.waterIce <1e-15   %resets STATUS back to zero
                snow.IA_PARENT.IA_PARENT_GROUND.IA_CHILD.STATUS = 0;
                snow = initialize_zero_snow(snow, snow.IA_PARENT.IA_PARENT_GROUND); %set all variables to zero
            end
            
            if snow.STATVAR.ice >= snow.PARA.swe_per_cell./2
                
                snow.IA_PARENT.IA_PARENT_GROUND.PREVIOUS.NEXT = snow; 
                snow.PREVIOUS = snow.IA_PARENT.IA_PARENT_GROUND.PREVIOUS;
                snow.NEXT = snow.IA_PARENT.IA_PARENT_GROUND;
                snow.IA_PARENT.IA_PARENT_GROUND.PREVIOUS = snow;
                snow.IA_PARENT.IA_PARENT_GROUND.IA_CHILD.STATUS = -1;
                snow.IA_PARENT.IA_PARENT_GROUND.IA_CHILD.FRACTIONAL_SNOW_COVER = 0;

                snow.IA_NEXT = get_IA_class(class(snow.NEXT), class(snow));
                snow.IA_PARENT.IA_PARENT_GROUND.IA_PREVIOUS = snow.IA_NEXT;
                snow.IA_NEXT.PREVIOUS = snow;
                snow.IA_NEXT.NEXT = snow.IA_PARENT.IA_PARENT_GROUND;
                
                %snow.IA_PARENT.IA_PARENT_GROUND.IA_CHILD.IA_CHILD_SNOW = [];
                snow.NEXT.IA_CHILD.IA_CHILD_SNOW = [];  %does not work yet to cut the connection between ground CHILD and snow
                snow.IA_PARENT = [];
            end
            % checks if snow CHILD needs to become full snow class and rearrange the stratigraphy
        end
        
        
        
        
        
        %non_mandatory functions
        
        function snow = conductivity(snow)   
            snow = conductivity@SNOW_base_class(snow);
        end
            
        function snow = get_T_water(snow)
            snow = get_T_water@SNOW_base_class(snow);
        end
        
        function snow = surface_energy_balance(snow, forcing)
            snow.STATVAR.Lout = (1-snow.PARA.epsilon) .* forcing.TEMP.Lin + snow.PARA.epsilon .* snow.CONST.sigma .* (snow.STATVAR.T(1)+ 273.15).^4;
            snow.STATVAR.Sout = snow.PARA.albedo .*  forcing.TEMP.Sin;
            snow.STATVAR.Qh = Q_h(snow, forcing);
            snow.STATVAR.Qe = Q_eq_potET(snow, forcing);
            
            snow.TEMP.F_ub = forcing.TEMP.Sin + forcing.TEMP.Lin - snow.STATVAR.Lout - snow.STATVAR.Sout - snow.STATVAR.Qh - snow.STATVAR.Qe;
        end
        
        
%         function snow = check_trigger(snow)
%             if size(snow.STATVAR.energy,1) ==1 && snow.STATVAR.ice < snow.PARA.swe_per_cell./2
%                 
%                 ground = snow.NEXT;
%                 
%                 ground.PREVIOUS = snow.PREVIOUS; %reassign ground
%                 ground.PREVIOUS.NEXT = ground;
%                 ground.IA_PREVIOUS=[];
%                 
%                 snow.NEXT =[]; %reassign snow -> cut all dependencies
%                 snow.PREVIOUS =[];
%                 snow.IA_NEXT =[];
%                 snow.IA_PREVIOUS =[];
%                 
%                 
%                 ground.IA_CHILD = IA_SNOW_GROUND();  %reinitialize interaction class
%                 ground.IA_CHILD.STATUS = 1; %snow initially active
%                 ground.IA_CHILD.IA_PARENT_GROUND = ground;  %attach snow and ground to interaction class
%                 ground.IA_CHILD.IA_CHILD_SNOW = snow;
%                 
%                 snow = ground; %assign snow pointer to ground to return to regular stratigraphy
%             end
%         end
        
        function snow = check_trigger(snow)
            if size(snow.STATVAR.energy,1) ==1 && snow.STATVAR.ice < snow.PARA.swe_per_cell/2
                
                ground = snow.NEXT;
                
                ground.PREVIOUS = snow.PREVIOUS; %reassign ground
                ground.PREVIOUS.NEXT = ground;
                ground.IA_PREVIOUS=[];
                
                snow.NEXT =[]; %reassign snow -> cut all dependencies
                snow.PREVIOUS =[];
                snow.IA_NEXT =[];
                snow.IA_PREVIOUS =[];
                
                %ground.IA_CHILD = IA_SNOW_GROUND_crocus();  %reinitialize interaction class
                ground.IA_CHILD.STATUS = 2; %snow initially active
                %ground.IA_CHILD.IA_PARENT_GROUND = ground;  %attach snow and ground to interaction class
                ground.IA_CHILD.IA_CHILD_SNOW = snow;
                ground.IA_CHILD.IA_CHILD_SNOW.IA_PARENT = ground.IA_CHILD;
                
                snow = ground; %assign snow pointer to ground to return to regular stratigraphy
            end
        end


    end
end