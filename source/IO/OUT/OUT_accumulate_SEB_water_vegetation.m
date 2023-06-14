%========================================================================
% CryoGrid OUT class defining storage format of the output
% OUT_accumulate_SEB_water_vegetation accumulates energy and water fluxes
% between out timesteps, for a user defined selection of classes. Also
% stores current layered vars (energy/water/temperature) of same classes.
% R. B. Zweigel, October 2022
%========================================================================


classdef OUT_accumulate_SEB_water_vegetation < matlab.mixin.Copyable
    
    
    properties
        out_index
        STATVAR
        TIMESTAMP
        HEIGHTS
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
        
    end
    
    
    methods
        
        %initialization
        
        function out = provide_PARA(out)
            
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
            out.PARA.classes = [];
        end
        
        
        function out = provide_CONST(out)
            
        end
        
        
        function out = provide_STATVAR(out)
            
            out.STATVAR.class   = [];
            % SEB fluxes
            out.STATVAR.Sin     = [];
            out.STATVAR.Lin     = [];
            out.STATVAR.Sout    = [];
            out.STATVAR.Lout    = [];
            out.STATVAR.Qe      = [];
            out.STATVAR.Qh      = [];
            % Water fluxes
            out.STATVAR.rainfall        = [];
            out.STATVAR.snowfall        = [];
            out.STATVAR.transp          = [];
            out.STATVAR.evap            = [];
            out.STATVAR.sublim          = [];
            out.STATVAR.surface_runoff  = [];
            out.STATVAR.subsurface_runoff = [];
            % Layer statvars
            out.STATVAR.T       = [];
            out.STATVAR.energy  = [];
            out.STATVAR.water   = [];
            out.STATVAR.ice     = [];
            out.STATVAR.waterIce = [];
        end
        
        
        function out = finalize_init(out, tile)
            % FINALIZE_SETUP  Performs all additional property
            %   initializations and modifications. Checks for some (but not
            %   all) data validity.
            
            %    ARGUMENTS:
            %    forcing:    instance of FORCING class
            forcing = tile.FORCING;
            
            % Set the next (first) output time. This is the next (first) time output
            % is collected (in memory) for later storage to disk.
            out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            
            % Set the next (first) save time. This is the next (first) time all the
            % collected output is saved to disk.
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval)
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
            
            out.TEMP = struct();
            out.TEMP.count = 0;
            out.TEMP.time = 0;
            
            % Energy
            out.TEMP.Sin    = zeros(size(out.PARA.classes));
            out.TEMP.Lin    = zeros(size(out.PARA.classes));
            out.TEMP.Sout   = zeros(size(out.PARA.classes));
            out.TEMP.Lout   = zeros(size(out.PARA.classes));
            out.TEMP.Qe     = zeros(size(out.PARA.classes));
            out.TEMP.Qh     = zeros(size(out.PARA.classes));
            % Water cycle
            out.TEMP.rainfall        = 0;
            out.TEMP.snowfall        = 0;
            out.TEMP.transp          = 0;
            out.TEMP.evap            = zeros(size(out.PARA.classes));
            out.TEMP.sublim          = zeros(size(out.PARA.classes));
            out.TEMP.surface_runoff  = 0;
            out.TEMP.subsurface_runoff = 0;
            
            heights = [];
            CURRENT = tile.TOP.NEXT;
            while ~isequal(CURRENT,tile.BOTTOM)
                for i = 1:length(out.PARA.classes)
                    className = class(CURRENT);
                    if strcmp(className,out.PARA.classes(i))
                        depths = [0; cumsum(CURRENT.STATVAR.layerThick)];
                        midpoints = CURRENT.STATVAR.upperPos - (depths(2:end)+depths(1:end-1))./2;
                        heights = [heights; midpoints];
                    end
                end
                
                CURRENT = CURRENT.NEXT;
            end
            
            if isfield(tile.STORE,'SNOW')
                zi = find(heights<=tile.PARA.altitude,1); %index of first ground cell
                heights = [heights(1:zi-1); [tile.PARA.altitude+1:-0.02:tile.PARA.altitude]'; heights(zi:end)]; % Snow out is hardcoded to 1m w. 2cm steps - Find smoother way later
            end
            
            out.HEIGHTS = heights;
        end
        
        %---------------time integration-------------
        
        %         function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
        
        function out = store_OUT(out, tile)
            
            
            t = tile.t;
            TOP = tile.TOP;
            BOTTOM = tile.BOTTOM;
            forcing = tile.FORCING;
            run_name = tile.PARA.run_name;
            result_path = tile.PARA.result_path;
            out_tag = out.PARA.tag;
            
            out.TEMP.count = out.TEMP.count + 1;
            out.TEMP.time = out.TEMP.time + tile.timestep;
            
            % ACCUMULATE VARIABLES
            % From forcing (but are downscaled???)
            
            CURRENT =TOP.NEXT;
            while ~isequal(CURRENT,BOTTOM)
                className = class(CURRENT);
                for i = 1:length(out.PARA.classes)
                    if strcmp(className(1:10),'VEGETATION') % Some vars only exist for vegetation
                        out.TEMP.transp = out.TEMP.transp + CURRENT.TEMP.transp .*tile.timestep;
                    end
                    if strcmp(className,out.PARA.classes(i))
                        % SEB
                        out.TEMP.Sin(i)  = out.TEMP.Sin(i) + CURRENT.STATVAR.Sin.*tile.timestep;
                        out.TEMP.Lin(i)  = out.TEMP.Lin(i) + CURRENT.STATVAR.Lin.*tile.timestep;
                        out.TEMP.Sout(i) = out.TEMP.Sout(i) + CURRENT.STATVAR.Sout.*tile.timestep;
                        out.TEMP.Lout(i) = out.TEMP.Lout(i) + CURRENT.STATVAR.Lout.*tile.timestep;
                        out.TEMP.Qe(i)   = out.TEMP.Qe(i) + CURRENT.STATVAR.Qe.*tile.timestep;
                        out.TEMP.Qh(i)   = out.TEMP.Qh(i) + CURRENT.STATVAR.Qh.*tile.timestep;
                        % Water
                        out.TEMP.evap(i)   = out.TEMP.evap(i) + CURRENT.TEMP.evap.*tile.timestep;
                        out.TEMP.sublim(i) = out.TEMP.sublim(i) + CURRENT.TEMP.sublim.*tile.timestep;
                    end
                end
                CURRENT = CURRENT.NEXT;
            end
            
            
            
            %           --- Write output ---
            if t>=out.OUTPUT_TIME
                %if id == 1
                avg_timestep = out.TEMP.time/out.TEMP.count;
                disp([datestr(t,'dd-mmm-yyyy HH:MM') ' Average timestep: ' num2str(avg_timestep)])
                out.TEMP.count = 0;
                out.TEMP.time = 0;
                %end
                %labBarrier
                
                out.TIMESTAMP=[out.TIMESTAMP t];
                
                out.TEMP.midpoints = [];
                out.TEMP.T = [];
                out.TEMP.waterIce = [];
                out.TEMP.energy = [];
                out.TEMP.water = [];
                out.TEMP.ice = [];
                
                CURRENT = TOP.NEXT;
                while ~isequal(CURRENT,BOTTOM)
                    for i = 1:length(out.PARA.classes)
                        className = class(CURRENT);
                        if strcmp(className,out.PARA.classes(i))
                            depths = [0; cumsum(CURRENT.STATVAR.layerThick)];
                            midpoints = CURRENT.STATVAR.upperPos - (depths(2:end)+depths(1:end-1))./2;
                            out.TEMP.midpoints = [out.TEMP.midpoints; midpoints; midpoints(end)-1e-4]; % add break (NaN) to separate classes
                            out.TEMP.T = [out.TEMP.T; CURRENT.STATVAR.T; NaN];
                            out.TEMP.waterIce = [out.TEMP.waterIce; CURRENT.STATVAR.waterIce; NaN];
                            out.TEMP.energy = [out.TEMP.energy; CURRENT.STATVAR.energy; NaN];
                            out.TEMP.water = [out.TEMP.water; CURRENT.STATVAR.water; NaN];
                            out.TEMP.ice = [out.TEMP.ice; CURRENT.STATVAR.ice; NaN];
                        end
                    end
                    CURRENT = CURRENT.NEXT;
                end
                
                % Energy
                out.STATVAR.Sin     = [out.STATVAR.Sin out.TEMP.Sin];
                out.STATVAR.Sout    = [out.STATVAR.Sout out.TEMP.Sout];
                out.STATVAR.Lin     = [out.STATVAR.Lin out.TEMP.Lin];
                out.STATVAR.Lout    = [out.STATVAR.Lout out.TEMP.Lout];
                out.STATVAR.Qe      = [out.STATVAR.Qe out.TEMP.Qe];
                out.STATVAR.Qh      = [out.STATVAR.Qh out.TEMP.Qh];
                % Water cycle
                out.STATVAR.transp  = [out.STATVAR.transp out.TEMP.transp];
                out.STATVAR.evap    = [out.STATVAR.evap out.TEMP.evap];
                out.STATVAR.sublim  = [out.STATVAR.sublim out.TEMP.sublim];
                % Layer statvars
                out.STATVAR.T       = [out.STATVAR.T interp1(out.TEMP.midpoints,out.TEMP.T,out.HEIGHTS)];
                out.STATVAR.waterIce= [out.STATVAR.waterIce interp1(out.TEMP.midpoints,out.TEMP.waterIce,out.HEIGHTS)];
                out.STATVAR.energy  = [out.STATVAR.energy interp1(out.TEMP.midpoints,out.TEMP.energy,out.HEIGHTS)];
                out.STATVAR.water   = [out.STATVAR.water interp1(out.TEMP.midpoints,out.TEMP.water,out.HEIGHTS)];
                out.STATVAR.ice     = [out.STATVAR.ice interp1(out.TEMP.midpoints,out.TEMP.ice,out.HEIGHTS)];
                
                out.TEMP.Sin    = out.TEMP.Sin.*0;
                out.TEMP.Lin    = out.TEMP.Lin.*0;
                out.TEMP.Sout   = out.TEMP.Sout.*0;
                out.TEMP.Lout   = out.TEMP.Lout.*0;
                out.TEMP.Qe     = out.TEMP.Qe.*0;
                out.TEMP.Qh     = out.TEMP.Qh.*0;
                out.TEMP.transp = out.TEMP.transp.*0;
                out.TEMP.evap   = out.TEMP.evap.*0;
                out.TEMP.sublim = out.TEMP.sublim.*0;
                
                % Set the next OUTPUT_TIME
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                
                
                if t>=out.SAVE_TIME
                    % It is time to save all the collected model output to disk
                    
                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                    else
                        save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                    end
                    
                    % Clear the out structure
                    out.TIMESTAMP=[];
                    out = provide_STATVAR(out);
                    if ~isnan(out.PARA.save_interval)
                        % If save_interval is defined, uptate SAVE_TIME for next save opertion
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        % If save_interval is not defined, we will save at the very end of the model run
                        % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
                        out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep); % OUTPUT_TIME is equal to SAVE_TIME when saving, need to update to avoid duplicating out of SAVE_TIME
                    end
                    
                end
            end
        end
        
        
    end
end