%%

classdef OUT_GEO4432
    properties
        TIMESTAMP
        RESULT
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
    end
    
    methods
        
        function xls_out = write_excel(out)
            xls_out = {'OUT','index',NaN,NaN;'OUT_GEO4432',1,NaN,NaN;'output_timestep',0.100000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'output_top',2,'[m]','heigh over ground included in output';'output_bottom',2,'[m]','depth under surface included in output';'output_spacing',0.0200000000000000,'[m]','vertical spacing og output';'OUT_END',NaN,NaN,NaN};
        end
        
        function out = provide_variables(out)
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.output_top = [];
            out.PARA.output_bottom = [];
            out.PARA.output_spacing = [];

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
            out.RESULT.grid = forcing.PARA.altitude + ...
                [-out.PARA.output_bottom:out.PARA.output_spacing:out.PARA.output_top]';

            out.TEMP.top_class = [];
            out.TEMP.time = 0;
            out.TEMP.Qh_acc = 0;
            out.TEMP.Qe_acc = 0;
            out.TEMP.Lin_acc = 0;
            out.TEMP.Lout_acc = 0;
            out.TEMP.Sin_acc = 0;
            out.TEMP.Sout_acc = 0;
            
            out.RESULT.Qh = [];
            out.RESULT.Qe = [];
            out.RESULT.Lin = [];
            out.RESULT.Lout = [];
            out.RESULT.Sin = [];
            out.RESULT.Sout = [];
            out.RESULT.T_surf = [];
            out.RESULT.d_snow = [];
            out.RESULT.swe = [];
            out.RESULT.T = [];
            
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval) ' 00:00:00'], 'dd.mm.yyyy HH:MM:SS'));
            end
        end
        
        function out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number, timestep, result_path)
            
            out.TEMP.Qh_acc = out.TEMP.Qh_acc - TOP_CLASS.STATVAR.Qh.*timestep;
            out.TEMP.Qe_acc = out.TEMP.Qe_acc - TOP_CLASS.STATVAR.Qe.*timestep;
            out.TEMP.Lout_acc = out.TEMP.Lout_acc - TOP_CLASS.STATVAR.Lout.*timestep;
            out.TEMP.Lin_acc = out.TEMP.Lin_acc + forcing.TEMP.Lin.*timestep;
            out.TEMP.Sout_acc = out.TEMP.Sout_acc - TOP_CLASS.STATVAR.Sout.*timestep;
            out.TEMP.Sin_acc = out.TEMP.Sin_acc + forcing.TEMP.Sin.*timestep;
            out.TEMP.time = out.TEMP.time + timestep;
            
            
            if t==out.OUTPUT_TIME
                disp([datestr(t)])
                out.TEMP.top_class = class(TOP_CLASS);
                out.TIMESTAMP=[out.TIMESTAMP t];
                
                if strcmp(out.TEMP.top_class(1:4),'SNOW')
                    out.RESULT.d_snow = [out.RESULT.d_snow sum(TOP_CLASS.STATVAR.layerThick)];
                    out.RESULT.swe = [out.RESULT.swe sum(TOP_CLASS.STATVAR.waterIce)];
                else
                    out.RESULT.d_snow = [out.RESULT.d_snow 0];
                    out.RESULT.swe = [out.RESULT.swe 0];
                end
                
                out.RESULT.T_surf = [out.RESULT.T_surf TOP_CLASS.STATVAR.T(1,1)];
                
                CURRENT = TOP_CLASS;
                layerThick=[];
                res = [];
                altitudeLowestCell = BOTTOM.PREVIOUS.STATVAR.lowerPos;
                while ~isequal(CURRENT, BOTTOM)
                    layerThick = [layerThick; CURRENT.STATVAR.layerThick];
                    res = [res; CURRENT.STATVAR.T];
                    CURRENT = CURRENT.NEXT;
                end
                % INCLUDE SOMETHING ABOUT SURFACE TEMP WHEN NO SNOW
                
                depths = [0; cumsum(layerThick)];
                depths = -(depths-depths(end,1));
                depths = (depths(1:end-1,1)+depths(2:end,1))./2 + altitudeLowestCell;
                
                %keyboard
                if max(depths(:)) < forcing.PARA.altitude
                    depths = [forcing.PARA.altitude; depths];
                    res = [TOP_CLASS.STATVAR.T(1,1); res];
                end
                
                result = interp1(depths, res, out.RESULT.grid);
                out.RESULT.T = [out.RESULT.T result ];
                    
                out.RESULT.Qh = [out.RESULT.Qh out.TEMP.Qh_acc./out.TEMP.time];
                out.RESULT.Qe = [out.RESULT.Qe out.TEMP.Qe_acc./out.TEMP.time];
                out.RESULT.Lout = [out.RESULT.Lout out.TEMP.Lout_acc./out.TEMP.time];
                out.RESULT.Lin = [out.RESULT.Lin out.TEMP.Lin_acc./out.TEMP.time];
                out.RESULT.Sout = [out.RESULT.Sout out.TEMP.Sout_acc./out.TEMP.time];
                out.RESULT.Sin = [out.RESULT.Sin out.TEMP.Sin_acc./out.TEMP.time];

                                
                out.TEMP.Qh_acc = 0;
                out.TEMP.Qe_acc = 0;
                out.TEMP.Lout_acc = 0;
                out.TEMP.Lin_acc = 0;
                out.TEMP.Sout_acc = 0;
                out.TEMP.Sin_acc = 0;
                out.TEMP.time = 0;

                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                
            
                if t >= out.SAVE_TIME 
                   
                   save([result_path run_number '/' run_number '_' datestr(t,'yyyy') '.mat'], 'out')
                   out.TIMESTAMP = [];
                   out.TEMP.top_class = [];
                   out.RESULT.Qh = [];
                   out.RESULT.Qe = [];
                   out.RESULT.Lin = [];
                   out.RESULT.Lout = [];
                   out.RESULT.Sin = [];
                   out.RESULT.Sout = [];
                   out.RESULT.T_surf = [];
                   out.RESULT.d_snow = [];
                   out.RESULT.swe = [];
                   out.RESULT.T = [];
                   out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                end
            end
        end
    end
    
end
