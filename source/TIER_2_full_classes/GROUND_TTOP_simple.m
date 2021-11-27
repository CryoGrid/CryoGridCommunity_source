%========================================================================
% CryoGrid GROUND class GROUND_TTOP_simple 
% S. Westermann, October 2020
%========================================================================

classdef GROUND_TTOP_simple < BASE
    
    methods
        
        
        function ground = provide_PARA(ground)

            ground.PARA.nt = []; %thawing n factor
            ground.PARA.nf = [];   %freezing n-factor
            ground.PARA.rk = [];    % rk
            ground.PARA.number_of_years = []; %in years, if empty, accumulate until end of time series
            ground.PARA.surface_T_variable = [];
            ground.PARA.timestep = [];
        end
        
        function ground = provide_STATVAR(ground)
            
            ground.STATVAR.TDD = [];
            ground.STATVAR.FDD = [];
            ground.STATVAR.MAGT = [];
            ground.STATVAR.MAGST = [];
            ground.STATVAR.MAAT = [];
            
        end
    
        function ground = provide_CONST(ground)
            ground.CONST.day_sec = [];

        end
        
        function ground = finalize_init(ground, tile)
            %ground.PARA.heatFlux_lb = tile.FORCING.PARA.heatFlux_lb;
            ground.TEMP.FDD_acc = 0;
            ground.TEMP.TDD_acc = 0;
            ground.TEMP.time_interval = 0;
            ground.TEMP.current_timestamp = tile.FORCING.start_time;
            ground.TEMP.next_timestamp = ground.TEMP.current_timestamp + ground.PARA.timestep;
            if isempty(ground.PARA.interval) || isnan(ground.PARA.interval)
                [~, month, day] = datevec(tile.FORCING.PARA.start_time);
                [year,~ ,~] = datevec(tile.FORCING.PARA.end_time);
            else
                [~, month, day] = datevec(tile.FORCING.PARA.start_time);
                [year,~ ,~] = datevec(tile.FORCING.PARA.start_time);
                year = year + ground.PARA.number_of_years;
            end
            ground.TEMP.next_TTOP_time = datenum(year, month, day);
            
        end
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile) 

        end
        
        function ground = get_boundary_condition_l(ground, tile)

        end
        
        function ground = get_derivatives_prognostic(ground, tile)

        end
        
        function timestep = get_timestep(ground, tile) 
           timestep = (min(ground.TEMP.next_timestamp, ground.TEMP.next_TTOP_time) - tile.t) .* ground.CONST.day_sec;
        end
        
        function ground = advance_prognostic(ground, tile) 
            timestep = tile.timestep;

            ground.TEMP.FDD_acc = ground.TEMP.FDD_acc + timestep ./ ground.CONST.day_sec .* tile.FORCING.TEMP.(ground.PARA.surface_T_variable) .* double(tile.FORCING.TEMP.(ground.PARA.surface_T_variable)<0);
            ground.TEMP.TDD_acc = ground.TEMP.TDD_acc + timestep ./ ground.CONST.day_sec .* tile.FORCING.TEMP.(ground.PARA.surface_T_variable) .* double(tile.FORCING.TEMP.(ground.PARA.surface_T_variable)>0);
            ground.TEMP.time_interval = ground.TEMP.time_interval + timestep ./ ground.CONST.day_sec;        
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)

        end
       
        function ground = compute_diagnostic(ground, tile)

            if tile.t + tile.timestep >= ground.TEMP.next_timestamp
                ground.TEMP.next_timestamp = ground.TEMP.next_timestamp + ground.PARA.timestep;
            end
            if tile.t + tile.timestep >= ground.TEMP.next_TTOP_time
                ground.STATVAR.MAAT = (ground.TEMP.TDD_acc + ground.TEMP.FDD_acc) ./ ground.TEMP.time_interval;

                ground.STATVAR.MAGST = (ground.PARA.nt .* ground.TEMP.TDD_acc + ground.PARA.nf .* ground.TEMP.FDD_acc) ./ ground.TEMP.time_interval;
                is_PF = ground.PARA.nt .* ground.PARA.rk .* ground.TEMP.TDD_acc + ground.PARA.nf  .* ground.TEMP.FDD_acc < 0;
                ground.STATVAR.MAGT = double(is_PF) .* (ground.PARA.nt .* ground.PARA.rk .* ground.TEMP.TDD_acc + ground.PARA.nf  .* ground.TEMP.FDD_acc) ./ ground.TEMP.time_interval +
                    double(~is_PF) .* (ground.PARA.nt .* ground.TEMP.TDD_acc + ground.PARA.nf ./ ground.PARA.rk  .* ground.TEMP.FDD_acc) ./ ground.TEMP.time_interval;

                ground.STATVAR.TDD = ground.TEMP.TDD_acc ./365.25;
                ground.STATVAR.FDD = ground.TEMP.FDD_acc ./365.25;

                ground.TEMP.FDD_acc = 0;
                ground.TEMP.TDD_acc = 0;
                ground.TEMP.time_interval = 0;
                
                [~, month, day] = datevec(tile.FORCING.PARA.start_time);
                [year,~ ,~] = datevec(tile.t);
                year = year + ground.PARA.number_of_years;
                ground.TEMP.next_TTOP_time = min(tile.FORCING.PARA.end_time, datenum(year, month, day));
            end
        end
        
        function ground = check_trigger(ground, tile)
           %do nothing 
        end
        

    end
    
end
