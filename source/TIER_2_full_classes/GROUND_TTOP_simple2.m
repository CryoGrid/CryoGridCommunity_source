%========================================================================
% CryoGrid GROUND class GROUND_TTOP_simple2
% GROUND_TTOP_simple2 computes MAAT, MAGST and MAAT based on the equilibrium
% TTOP approach.
% S. Westermann, October 2020
%========================================================================

classdef GROUND_TTOP_simple2 < BASE
    
    methods
        
        
        function ground = provide_PARA(ground)

            ground.PARA.nt = []; %thawing n factor
            ground.PARA.nf = [];   %freezing n-factor
            ground.PARA.rk = [];    % rk
            ground.PARA.number_of_years = []; %in years, if empty, accumulate until end of time series
            ground.PARA.start_date = [];
            ground.PARA.surface_T_variable = [];
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.MAGT = [];
            ground.STATVAR.MAGST = [];
            ground.STATVAR.MAAT = [];
            
        end
    
        function ground = provide_CONST(ground)
            ground.CONST.day_sec = [];

        end
        
        function ground = finalize_init(ground, tile)
            %ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            timestamp = tile.FORCING.DATA.timeForcing;
            surface_T = tile.FORCING.DATA.(ground.PARA.surface_T_variable);
            ground.STATVAR.layerThick = 0;
            
            if tile.FORCING.PARA.start_time<= datenum([datestr(tile.FORCING.PARA.start_time, 'yyyy') ground.PARA.start_date], 'yyyymmdd')
                start_time = datenum([datestr(tile.FORCING.PARA.start_time, 'yyyy') ground.PARA.start_date], 'yyyymmdd');
            else
                start_time = datenum([num2str(year(tile.FORCING.PARA.start_time)+1) ground.PARA.start_date], 'yyyymmdd');
            end
                
            if tile.FORCING.PARA.end_time>= datenum([datestr(tile.FORCING.PARA.end_time, 'yyyy') ground.PARA.start_date], 'yyyymmdd')
                end_time = datenum([datestr(tile.FORCING.PARA.end_time, 'yyyy') ground.PARA.start_date], 'yyyymmdd');
            else
                end_time = datenum([num2str(year(tile.FORCING.PARA.end_time)-1) ground.PARA.start_date], 'yyyymmdd');
            end

            ground.TEMP.end_timestamp = [];
            ground.TEMP.start_timestamp = start_time;
            while datenum([num2str(year(start_time) + ground.PARA.number_of_years) ground.PARA.start_date], 'yyyymmdd') < end_time
                start_time = datenum([num2str(year(start_time) + ground.PARA.number_of_years) ground.PARA.start_date], 'yyyymmdd');
                ground.TEMP.end_timestamp = [ground.TEMP.end_timestamp; start_time];
                ground.TEMP.start_timestamp = [ground.TEMP.start_timestamp; start_time];
            end
            ground.TEMP.end_timestamp = [ground.TEMP.end_timestamp; end_time];
            
            ground.TEMP.FDD = [];
            ground.TEMP.TDD = [];
            number_of_values = [];
            
            for i=1:size(ground.TEMP.start_timestamp, 1)
                range = find(timestamp(:,1)>=ground.TEMP.start_timestamp(i,1) & timestamp(:,1)<ground.TEMP.end_timestamp(i,1));
                temperature = surface_T(range,1);
                ground.TEMP.FDD = [ground.TEMP.FDD; sum(temperature .* double(temperature <0)) ];
                ground.TEMP.TDD = [ground.TEMP.TDD; sum(temperature .* double(temperature >0)) ];
                number_of_values = [number_of_values; length(range)];
            end
            ground.TEMP.end_timestamp(end) = tile.FORCING.PARA.end_time;
            
            ground.TEMP.MAAT = (ground.TEMP.FDD + ground.TEMP.TDD) ./ number_of_values;
            ground.TEMP.MAGST = (ground.PARA.nf .* ground.TEMP.FDD + ground.PARA.nt .* ground.TEMP.TDD) ./ number_of_values;           
            is_PF = ground.PARA.nt .* ground.PARA.rk .* ground.TEMP.TDD + ground.PARA.nf  .* ground.TEMP.FDD < 0;
            ground.TEMP.MAGT = double(is_PF) .* (ground.PARA.nt .* ground.PARA.rk .* ground.TEMP.TDD + ground.PARA.nf  .* ground.TEMP.FDD) ./ number_of_values + ...
                    double(~is_PF) .* (ground.PARA.nt .* ground.TEMP.TDD + ground.PARA.nf ./ ground.PARA.rk  .* ground.TEMP.FDD) ./ number_of_values;
            ground.TEMP.interval = 1;
            
            ground.STATVAR.MAAT = ground.TEMP.MAAT(ground.TEMP.interval,1);
            ground.STATVAR.MAGST = ground.TEMP.MAGST(ground.TEMP.interval,1);
            ground.STATVAR.MAGT = ground.TEMP.MAGT(ground.TEMP.interval,1);
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile) 
            
        end
        
        function ground = get_boundary_condition_l(ground, tile)

        end
        
        function ground = get_derivatives_prognostic(ground, tile)

        end
        
        function timestep = get_timestep(ground, tile) 
           timestep = (ground.TEMP.end_timestamp(ground.TEMP.interval,1) - tile.t) .* ground.CONST.day_sec;
        end
        
        function ground = advance_prognostic(ground, tile) 
    
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)

        end
       
        function ground = compute_diagnostic(ground, tile)
            if tile.t >= ground.TEMP.end_timestamp(ground.TEMP.interval,1)
                ground.TEMP.interval = ground.TEMP.interval + 1;
            end
            ground.STATVAR.MAAT = ground.TEMP.MAAT(ground.TEMP.interval,1);
            ground.STATVAR.MAGST = ground.TEMP.MAGST(ground.TEMP.interval,1);
            ground.STATVAR.MAGT = ground.TEMP.MAGT(ground.TEMP.interval,1);
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end
        
        
        %-------------param file generation-----
        function ground = param_file_info(ground)
            ground = param_file_info@BASE(ground);
            
            ground.PARA.class_category = 'GROUND';
            
            ground.PARA.STATVAR = {''};
            
            ground.PARA.default_value.nt = {1};
            ground.PARA.comment.nt = {'thawing n-factor'};
            
            ground.PARA.default_value.nf = {0.5};
            ground.PARA.comment.nf = {'freezing n-factor'};
            
            ground.PARA.default_value.rk = {0.8};
            ground.PARA.comment.rk = {'ratio thawed to frozen thermal conductivity'};
            
            ground.PARA.default_value.number_of_years = {10};
            ground.PARA.comment.number_of_years = {'number of years to calculate constant TTOP value' };
            
            ground.PARA.default_value.surface_T_variable = {'Tair'};
            ground.PARA.comment.surface_T_variable = {'variable in FORCING to be used as surfac temperature'};
            
            ground.PARA.default_value.start_date = {'0101'};
            ground.PARA.comment.start_date = {'date of the year mmdd to start TTOP intervals'};
            
        end

    end
    
end
