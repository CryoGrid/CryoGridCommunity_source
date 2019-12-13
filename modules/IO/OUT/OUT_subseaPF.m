%> OUT_SUBSEAPF Manage and save output data for subsea permafrost
%> manages a struct with forcing data, info on the performance of the
%> current run and stores the results at certain times
%> saves the accumulated results at the final time or if the run breaks
classdef OUT_subseaPF
    properties
        STRATIGRAPHY
        TIMESTAMP
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        
        FORCING
        BREAK
        RUNINFO %save: anzahl echte zeitschritte, rechenzeit (gesamt)
    end
    
    
    methods
              
        function out = provide_variables(out)
            out.PARA.output_timestep = [];
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
           out.PARA.output_timestep = out.PARA.output_timestep*365.24;
        end
        
        function out = complete_init_out(out, forcing)
            out.OUTPUT_TIME = (forcing.PARA.startForcing)*365.24;%start saving on the first time step + out.PARA.output_timestep;
            out.SAVE_TIME = (forcing.PARA.endForcing)*365.24;
            
            %save forcing data
            out.FORCING.airTemp = forcing.DATA.airTemp;
            out.FORCING.seaLevel = forcing.DATA.seaLevel;
            out.FORCING.glacialCover = forcing.DATA.glacialCover;
            out.FORCING.timeForcing = forcing.DATA.timeForcing;
            out.FORCING.TForcing = [];
            out.FORCING.saltConcForcing = [];
            
            %initialise Runtime and timesteps
            out.RUNINFO.starttime = tic;
            out.RUNINFO.timesteps = 0;
            out.RUNINFO.dt_min = 1000;
            out.RUNINFO.dt_max = 0;
            
            out.PARA.lastDisp = out.OUTPUT_TIME -100;
            out.BREAK = 0;
            
            figure(1)
            clf
        end
        
        function out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number)
            
         
            if isnan(mean(TOP_CLASS.STATVAR.T))
                fprintf('Time is %.0f - temperature is NAN - terminating! \n', t/365.25)
                out.BREAK = 1;
            end
            if isfield(TOP_CLASS.STATVAR, 'saltConc') && isnan(mean(TOP_CLASS.STATVAR.saltConc))
                fprintf('Time is %.0f - salt concentration is NAN - terminating! \n', t/365.25)
                out.BREAK = 1;
            end
            
            %display current state every now and then
            dispInterval = 100*365.25;
            if abs(out.PARA.lastDisp - t) > dispInterval
                out.PARA.lastDisp = t;
                fprintf('Time is %.0f years bfi \n', t/365.25)
                            
                if 0%~isempty(out.STRATIGRAPHY) % && stratigraphy{1,1}.TEMP.saltFlux_ub ~= 0 %something
                    stratigraphy = out.STRATIGRAPHY{1,end};

                    figure(1)
                    subplot(3,1,1)
                    yyaxis left
                    plot(t/365.25, stratigraphy{1,1}.TEMP.T_ub, 'o')
                    hold on
                    if isfield(stratigraphy{1,1}.TEMP, 'saltConc_ub')
                        yyaxis right
                        plot(t/365.25, stratigraphy{1,1}.TEMP.saltConc_ub, 'x')
                    end
                    drawnow

                    subplot(3,1,[2:3])
                    CURRENT = TOP_CLASS;
                    while ~isequal(CURRENT, BOTTOM)
                        if isprop(CURRENT, 'IA_CHILD') && ~isempty(CURRENT.IA_CHILD) && CURRENT.IA_CHILD.STATUS ==1
                            plot(CURRENT.IA_CHILD.IA_CHILD.STATVAR.T, CURRENT.IA_CHILD.IA_CHILD.STATVAR.upperPos - cumsum(CURRENT.IA_CHILD.IA_CHILD.STATVAR.layerThick) - 0.5*CURRENT.IA_CHILD.IA_CHILD.STATVAR.layerThick)
                        end
                    
                        plot(CURRENT.STATVAR.T, CURRENT.STATVAR.upperPos - cumsum(CURRENT.STATVAR.layerThick) - 0.5*CURRENT.STATVAR.layerThick)
                        hold on
    %                     fprintf('Upper Position: %3.1f \nTemperature: %2.1f\n', CURRENT.STATVAR.upperPos, CURRENT.STATVAR.T(1))
    %                     fprintf('Lower Position: %3.1f \nTemperature: %2.1f\n', CURRENT.STATVAR.lowerPos, CURRENT.STATVAR.T(end))
                        CURRENT = CURRENT.NEXT;
                    end
                    ylim(stratigraphy{1,1}.STATVAR.upperPos + [-2, 0])
                    drawnow
                end
            end
            
            %update runinfo every step
            out.RUNINFO.timesteps = out.RUNINFO.timesteps + 1;
            %out.RUNINFO.dt_min = min(out.RUNINFO.dt_min, run_info.current_timestep) ;
            %out.RUNINFO.dt_max = max(out.RUNINFO.dt_max, run_info.current_timestep);
            
            
            if t==out.OUTPUT_TIME || out.BREAK == 1
                %save forcing
                out.FORCING.TForcing = [out.FORCING.TForcing, forcing.TEMP.TForcing];
                out.FORCING.saltConcForcing = [out.FORCING.saltConcForcing, forcing.TEMP.saltConcForcing];
                            
                %save sediment temperature
                out.TIMESTAMP=[out.TIMESTAMP t];
                CURRENT =TOP_CLASS;
                if isprop(CURRENT, 'IA_CHILD') && ~isempty(CURRENT.IA_CHILD)
                    out.MISC=[out.MISC [CURRENT.IA_CHILD.IA_CHILD.STATVAR.T(1,1); CURRENT.IA_CHILD.IA_CHILD.STATVAR.layerThick(1,1)]]; 
                else
                    out.MISC=[out.MISC [NaN; NaN]];
                end
                result={}; 
                                
                while ~isequal(CURRENT, BOTTOM)
                    if isprop(CURRENT, 'IA_CHILD') && ~isempty(CURRENT.IA_CHILD) && CURRENT.IA_CHILD.STATUS ==1
                        res=copy(CURRENT.IA_CHILD.IA_CHILD);
                        res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  %cut all dependencies
                        result=[result; {res}];
                    end
                    res = copy(CURRENT);
                    res.NEXT =[]; res.PREVIOUS=[]; res.IA_NEXT=[]; res.IA_NEXT=[];  %cut all dependencies
                    if isprop(res, 'IA_CHILD')
                        res.IA_CHILD =[];
                    end
                    result=[result; {res}];
                    CURRENT = CURRENT.NEXT;
                end
                out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = result;
                
                        
                %make sure to hit savetime at the end
                out.OUTPUT_TIME = min(out.OUTPUT_TIME + out.PARA.output_timestep, out.SAVE_TIME);            
                
                if t==out.SAVE_TIME || out.BREAK == 1
                    
                    out.RUNINFO.endtime = toc(out.RUNINFO.starttime);
                    save(['results' filesep run_number filesep 'Results_' run_number '.mat'], 'out', '-v7.3')
                    out.BREAK = 1;

                end
            end
        end
        
    end
end