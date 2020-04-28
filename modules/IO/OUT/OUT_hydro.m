% OUT class which outputs data ready to be plotted
% R. Zweigel, October 2019
classdef OUT_hydro
    properties
        TIMESTAMP
        RESULT
        META
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
    end
    
    methods
        
        function xls_out = write_excel(out)
            xls_out = {'OUT','index',NaN,NaN;'OUT_parallel',1,NaN,NaN;'output_timestep',0.125000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'height',1,'[m]','Height above ground included';'depth',1,'[m]','depth below ground included';'ground_spacing',0.0250000000000000,'[m]','interval below ground';'snow_spacing',0.0100000000000000,'[m]','interval in snow';'status_seb',1,NaN,'1 = on, 0 = off';'status_snow',1,NaN,'1 = on, 0 = off';'OUT_END',NaN,NaN,NaN};
        end
        
        function out = provide_variables(out)
            out.PARA.output_timestep    = [];
            out.PARA.save_date          = [];
            out.PARA.save_interval      = [];
            out.PARA.height             = [];
            out.PARA.depth              = [];
            out.PARA.ground_spacing     = [];
            out.PARA.gst_depth          = [];
            out.PARA.snow_spacing       = [];
            out.PARA.status_hydro       = [];
            out.PARA.status_snow        = [];
            out.META.altitude           = [];
            out.META.longitude          = [];
            out.META.latitude           = [];
            out.META.forcing_name       = [];
        end
        
        function out = initalize_from_file(out, section)
            variables = fieldnames(out.PARA);
            for i=1:size(variables,1)
                for j=1:size(section,1)
                    if strcmp(variables{i,1}, section{j,1})
                        out.PARA.(variables{i,1}) = section{j,2};
                    end
                end
            end
        end
        
        function out = complete_init_out(out, forcing)
            out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            
            if out.PARA.status_snow >= 1
                out.RESULT.grid_snow = out.PARA.height:-out.PARA.snow_spacing:0;
                out.RESULT.grid_snow = double(out.RESULT.grid_snow + forcing.PARA.altitude);
            end
            
            out.RESULT.grid = [out.PARA.height:-out.PARA.snow_spacing:0,...
                -out.PARA.ground_spacing:-out.PARA.ground_spacing:-out.PARA.depth];
            out.RESULT.grid = double(out.RESULT.grid + forcing.PARA.altitude);
            
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval)
                out.SAVE_TIME = floor(forcing.PARA.end_time);
            else
                out.SAVE_TIME = min(floor(forcing.PARA.end_time),  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval) ' 00:00:00'], 'dd.mm.yyyy HH:MM:SS'));
            end
            
            % Metadata of current run
            out.META.altitude           = forcing.PARA.altitude;
            out.META.longitude          = forcing.PARA.longitude;
            out.META.latitude           = forcing.PARA.latitude;
            out.META.forcing_name       = forcing.PARA.filename;
            out.META.rain_fraction      = forcing.PARA.rain_fraction;
            out.META.snow_fraction      = forcing.PARA.snow_fraction;
            
            % Basic variables
            out.RESULT.T            = [];
            out.RESULT.water        = [];
            out.RESULT.ice          = [];
            out.RESULT.FDD          = 0;
            out.RESULT.TDD          = 0;
            
            
            if out.PARA.status_hydro == 1
                out.TEMP.top_class  = [];
                out.TEMP.time       = 0;
                out.TEMP.ET_acc     = 0;
                out.TEMP.precp_acc  = 0;
                out.TEMP.lateral_water = 0;
                
                out.RESULT.runoff   = [];
                out.RESULT.precip   = [];
                out.RESULT.ET       = [];
                out.RESULT.lateral  = [];
                out.RESULT.storage  = [];
            end
            if out.PARA.status_snow == 1
                out.RESULT.swe      = [];
                out.RESULT.d_snow   = [];
                out.RESULT.saturation = [];
                out.RESULT.lwc      = [];
                out.RESULT.density  = [];
            end
            if out.PARA.status_snow == 2
                out.RESULT.swe      = [];
                out.RESULT.d_snow   = [];
                out.RESULT.d        = [];
                out.RESULT.s        = [];
                out.RESULT.gs       = [];
                out.RESULT.saturation = [];
                out.RESULT.lwc      = [];
                out.RESULT.density  = [];
                out.RESULT.albedo       = [];
            end
            
        end
        
        function out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number, timestep, result_path, lateral)
            
            if out.PARA.status_hydro == 1
                % Storage initial state
                if isempty(out.TEMP.top_class)
                    CURRENT = TOP_CLASS;
                    out.RESULT.storage = 0;
                    while ~isequal(CURRENT, BOTTOM)
                        out.RESULT.storage = out.RESULT.storage + sum(CURRENT.STATVAR.waterIce);
                        CURRENT = CURRENT.NEXT;
                    end
                end
                
                out.TEMP.top_class = class(TOP_CLASS);
                
                % EvapoTranspiration
                if strcmp(out.TEMP.top_class(1:6),'GROUND')
                    L = 1e3.*(2500.8 - 2.36.*TOP_CLASS.STATVAR.T(1)) .*1000;
                    out.TEMP.ET_acc     = out.TEMP.ET_acc + -TOP_CLASS.STATVAR.Qe ./ L*timestep;
                    
                elseif strcmp(out.TEMP.top_class(1:4),'SNOW')
                    out.TEMP.ET_acc     = out.TEMP.ET_acc + -TOP_CLASS.STATVAR.Qe ./ TOP_CLASS.CONST.L_s*timestep;
                end
                
                out.TEMP.precp_acc  = out.TEMP.precp_acc + (forcing.TEMP.rainfall + forcing.TEMP.snowfall)./1000 ./ (24.*3600)*timestep;
                
                if numlabs > 1 % lateral fluxes can occur
                    out.TEMP.lateral_water = out.TEMP.lateral_water + TOP_CLASS.TEMP.lateral_water;
                end
                out.TEMP.time = out.TEMP.time + timestep;
                
            end
            
            if t==out.OUTPUT_TIME
                if exist('lateral') && labindex == 1
                    disp([datestr(t,'dd-mmm-yyyy HH:MM:SS') ' lateral status; snow: ' num2str(lateral.STATUS.snow) ' water: ' num2str(lateral.STATUS.water) ' waterflux: ' num2str(out.TEMP.lateral_water)])
                elseif ~exist('lateral')
                    disp(datestr(t,'dd-mmm-yyyy HH:MM:SS'))
                end
                
                out.TIMESTAMP=[out.TIMESTAMP t];
                
                CURRENT = TOP_CLASS;
                layerThick  = [];
                T           = [];
                water       = [];
                ice         = [];
                while ~isequal(CURRENT, BOTTOM)
                    layerThick = [layerThick; CURRENT.STATVAR.layerThick];
                    T = [T; CURRENT.STATVAR.T];
                    water = [water; CURRENT.STATVAR.water];
                    ice = [ice; CURRENT.STATVAR.ice];
                    CURRENT = CURRENT.NEXT;
                end
                depths = [0; cumsum(layerThick)];
                depths = -(depths-depths(end,1));
                depths = (depths(1:end-1,1)+depths(2:end,1))./2 + BOTTOM.PREVIOUS.STATVAR.lowerPos;
                depths = [TOP_CLASS.STATVAR.upperPos; depths];
                
                water   = water ./ layerThick;
                ice     = ice ./ layerThick;
                
                T       = [T(1); T];
                water   = [water(1); water];
                ice     = [ice(1); ice];
                
                if  out.PARA.status_snow >= 1
                    if strcmp(out.TEMP.top_class(1:4),'SNOW')
                        out.RESULT.swe      = [out.RESULT.swe sum(TOP_CLASS.STATVAR.waterIce)];
                        out.RESULT.d_snow   = [out.RESULT.d_snow sum(TOP_CLASS.STATVAR.layerThick)];
                        out.RESULT.lwc      = [out.RESULT.lwc sum(TOP_CLASS.STATVAR.water)];
                        snowdepths  = [0; cumsum(TOP_CLASS.STATVAR.layerThick)];
                        snowdepths  = -(snowdepths - snowdepths(end));
                        snowdepths  = (snowdepths(1:end-1) + snowdepths(2:end))./2 + TOP_CLASS.STATVAR.lowerPos;
                        snowdepths  = [TOP_CLASS.STATVAR.upperPos; snowdepths; TOP_CLASS.STATVAR.lowerPos];
                        if out.PARA.status_snow == 2
                            d  = TOP_CLASS.STATVAR.d;
                            s  = TOP_CLASS.STATVAR.s;
                            gs = TOP_CLASS.STATVAR.gs;
                            d  = [d(1); d; d(end)];
                            s  = [s(1); s; s(end)];
                            gs = [gs(1); gs; gs(end)];
                            albedo = TOP_CLASS.TEMP.albedo;
                        end
                        porespace   = TOP_CLASS.STATVAR.layerThick - TOP_CLASS.STATVAR.ice;
                        saturation  = TOP_CLASS.STATVAR.water./porespace;
                        saturation(porespace == 0) = 0;
                        free_water  = max(0, TOP_CLASS.STATVAR.water - (TOP_CLASS.STATVAR.layerThick - TOP_CLASS.STATVAR.ice).*TOP_CLASS.PARA.field_capacity);
                        density     = (TOP_CLASS.STATVAR.ice.*917 + (TOP_CLASS.STATVAR.water-free_water).*1000)./TOP_CLASS.STATVAR.layerThick;
                        saturation  = [saturation(1); saturation; saturation(end)];
                        density     = [density(1); density; density(end)];
                        
                        if sum(porespace == 0 & TOP_CLASS.STATVAR.water ~= 0) > 0
                            dfe4w
                        end
                        
                    elseif TOP_CLASS.IA_CHILD.STATUS == 1 && TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.layerThick > eps(TOP_CLASS.STATVAR.upperPos)
                        out.RESULT.d_snow   = [out.RESULT.d_snow TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.layerThick];
                        out.RESULT.swe      = [out.RESULT.swe TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.waterIce];
                        out.RESULT.lwc      = [out.RESULT.lwc TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.water];
                        snowdepths  = [TOP_CLASS.STATVAR.upperPos + TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.layerThick; TOP_CLASS.STATVAR.upperPos];
                        if out.PARA.status_snow == 2
                            d   = repmat(TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.d,2,1);
                            s   = repmat(TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.s,2,1);
                            gs  = repmat(TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.gs,2,1);
                            albedo = TOP_CLASS.TEMP.albedo;
                        end
                        porespace   = (TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.layerThick - TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.ice);
                        saturation  = repmat(TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.water./porespace,2,1);
                        saturation(porespace == 0) = 0;
                        density     = repmat(TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.target_density,2,1); % RBZ: Pulling a fas one - not consistent method
                    else
                        out.RESULT.d_snow   = [out.RESULT.d_snow 0];
                        out.RESULT.swe      = [out.RESULT.swe 0];
                        out.RESULT.lwc      = [out.RESULT.lwc 0];
                        snowdepths  = [TOP_CLASS.STATVAR.upperPos+.01; TOP_CLASS.STATVAR.upperPos];
                        if out.PARA.status_snow == 2
                            d   = [NaN; NaN];
                            s   = [NaN; NaN];
                            gs  = [NaN; NaN];
                            albedo = TOP_CLASS.PARA.albedo;
                        end
                        density     = [NaN; NaN];
                        saturation  = [NaN; NaN];
                    end
                end
                
                out.RESULT.T        = [out.RESULT.T interp1(depths, T, out.RESULT.grid)' ];
                out.RESULT.water    = [out.RESULT.water interp1(depths, water, out.RESULT.grid)'];
                out.RESULT.ice      = [out.RESULT.ice interp1(depths, ice, out.RESULT.grid)' ];
                
                % TDD/FDD calculation
                GST                 = interp1(depths, T, out.META.altitude - out.PARA.gst_depth);
                if GST > 0
                    out.RESULT.TDD      = out.RESULT.TDD + GST * out.PARA.output_timestep;
                elseif GST < 0
                    out.RESULT.FDD      = out.RESULT.FDD + GST * out.PARA.output_timestep;
                end
                
                if  out.PARA.status_snow >= 1
                    if out.PARA.status_snow == 2
                        out.RESULT.d        = [out.RESULT.d interp1(snowdepths, d, out.RESULT.grid_snow)' ];
                        out.RESULT.s        = [out.RESULT.s interp1(snowdepths, s, out.RESULT.grid_snow)' ];
                        out.RESULT.gs       = [out.RESULT.gs interp1(snowdepths, gs, out.RESULT.grid_snow)' ];
                        out.RESULT.albedo   = [out.RESULT.albedo albedo];
                    end
                    out.RESULT.saturation   = [out.RESULT.saturation interp1(snowdepths, saturation, out.RESULT.grid_snow)' ];
                    out.RESULT.density      = [out.RESULT.density interp1(snowdepths, density, out.RESULT.grid_snow)' ];
                end
                
                if out.PARA.status_hydro == 1
                    
                    % Add the water in the CHILD
                    if strcmp(out.TEMP.top_class(1:6),'GROUND')
                        storage = TOP_CLASS.IA_CHILD.IA_CHILD_SNOW.STATVAR.waterIce;
                    else
                        storage = 0;
                    end
                    % Add water from classes
                    CURRENT = TOP_CLASS;
                    while ~isequal(CURRENT, BOTTOM)
                        storage = storage + sum(CURRENT.STATVAR.waterIce);
                    	CURRENT = CURRENT.NEXT;
                    end
                    
                    out.RESULT.precip   = [out.RESULT.precip out.TEMP.precp_acc./out.TEMP.time];
                    out.RESULT.ET       = [out.RESULT.ET out.TEMP.ET_acc./out.TEMP.time];
                    out.RESULT.storage  = [out.RESULT.storage storage];
                    out.RESULT.lateral  = [out.RESULT.lateral out.TEMP.lateral_water./out.TEMP.time]; 
                    out.RESULT.runoff   = [out.RESULT.runoff out.RESULT.precip(end)+ out.RESULT.ET(end)...
                        + out.RESULT.lateral(end) + (out.RESULT.storage(end-1)-out.RESULT.storage(end))./out.TEMP.time];

                    out.TEMP.precp_acc      = 0;
                    out.TEMP.ET_acc         = 0;
                    out.TEMP.lateral_water  = 0;
                    out.TEMP.time           = 0;
                end
                
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                
                if t >= out.SAVE_TIME
                    
                    save([result_path run_number '/' run_number num2str(labindex) '_' datestr(t,'yyyy') '.mat'], 'out')
                    out.TIMESTAMP       = [];
                    out.TEMP.top_class  = [];
                    out.RESULT.T        = [];
                    out.RESULT.water    = [];
                    out.RESULT.ice      = [];
                    out.RESULT.FDD      = 0;
                    out.RESULT.TDD      = 0;
                    
                    if out.PARA.status_hydro == 1
                        out.RESULT.precip   = [];
                        out.RESULT.ET       = [];
                        out.RESULT.runoff   = [];
                        out.RESULT.lateral  = [];
                        out.RESULT.storage  = out.RESULT.storage(end);

                    end
                    if out.PARA.status_snow >= 1
                        out.RESULT.d_snow   = [];
                        out.RESULT.swe      = [];
                        out.RESULT.lwc      = [];
                        if out.PARA.status_snow == 2
                            out.RESULT.d        = [];
                            out.RESULT.s        = [];
                            out.RESULT.gs       = [];
                            out.RESULT.albedo   = [];
                        end
                        out.RESULT.saturation = [];
                        out.RESULT.density  = [];
                    end
                    
                    out.SAVE_TIME = min(floor(forcing.PARA.end_time),  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end
        end
        
    end
end
