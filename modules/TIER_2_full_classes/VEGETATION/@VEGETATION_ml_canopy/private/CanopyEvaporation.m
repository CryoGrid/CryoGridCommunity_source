function [vegetation] = CanopyEvaporation (vegetation, p) %, ic, il

% %DESCRIPTION:
% Update canopy intercepted water for evaporation and dew

% %USES:
% dtime_sub = vegetation.physcon.dtime_sub;                  %clm_varctl
isun = vegetation.params.sun;   %Sunlit leaf
isha = vegetation.params.sha;   %Shaded leaf
% mmh2o = vegetation.physcon.mmh2o;                       %use clm_varcon

% p = vegetation.mlcanopyinst.p;
% ncan = vegetation.mlcanopyinst.ncan;
% vegetation.canopy.dpai = vegetation.mlcanopyinst.dpai;

%--------------------------------------------------------------------
        % 1.2 CanopyEvaporation
%Input 
%trleaf = vegetation.mlcanopyinst.trleaf;                          % mlcanopy_inst%trleaf              % Leaf transpiration flux (mol H2O/m2 leaf/s)
%vegetation.mlcanopyinst.evleaf = vegetation.mlcanopyinst.evleaf;                          % mlcanopy_inst%evleaf              % Leaf evaporation flux (mol H2O/m2 leaf/s)
%vegetation.flux.fracsun = vegetation.flux.fracsun;                        % mlcanopy_inst%fracsun             % Sunlit fraction of canopy layer
%vegetation.flux.fracsha = vegetation.flux.fracsha;                        % mlcanopy_inst%vegetation.flux.fracsha             % Shaded fraction of canopy layer

%Input/Output
%vegetation.mlcanopyinst.h2ocan = vegetation.mlcanopyinst.h2ocan;                        % mlcanopy_inst%h2ocan              % Canopy layer intercepted water (kg H2O/m2)
%--------------------------------------------------------------------

    dtime = vegetation.params.dtime_sub;
    
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
        
        if (vegetation.canopy.dpai(p,ic) > 0)
            
            % Add dew, from both evaporation and transpiration
            
            dew = (vegetation.mlcanopyinst.evleaf(p,ic,isun) + vegetation.mlcanopyinst.trleaf(p,ic,isun)) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic) * vegetation.physcon.mmh2o * dtime;
            if (dew < 0.)
                vegetation.mlcanopyinst.h2ocan(p,ic) = vegetation.mlcanopyinst.h2ocan(p,ic) - dew;
            end
            
            dew = (vegetation.mlcanopyinst.evleaf(p,ic,isha) + vegetation.mlcanopyinst.trleaf(p,ic,isha)) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic) * vegetation.physcon.mmh2o * dtime;
            if (dew < 0.)
                vegetation.mlcanopyinst.h2ocan(p,ic) = vegetation.mlcanopyinst.h2ocan(p,ic) - dew;
            end
            
            % Evaporate intercepted water
            
            if (vegetation.mlcanopyinst.evleaf(p,ic,isun) > 0)
                vegetation.mlcanopyinst.h2ocan(p,ic) = vegetation.mlcanopyinst.h2ocan(p,ic) - vegetation.mlcanopyinst.evleaf(p,ic,isun) * vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic) * vegetation.physcon.mmh2o * dtime;
            end
            
            if (vegetation.mlcanopyinst.evleaf(p,ic,isha) > 0)
                vegetation.mlcanopyinst.h2ocan(p,ic) = vegetation.mlcanopyinst.h2ocan(p,ic) - vegetation.mlcanopyinst.evleaf(p,ic,isha) * vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic) * vegetation.physcon.mmh2o * dtime;
            end
            
            % for not allow vegetation.mlcanopyinst.h2ocan to go negative
            % vegetation.mlcanopyinst.h2ocan(p,ic) = max (0._r8, vegetation.mlcanopyinst.h2ocan(p,ic))
            
            % Exclude evaporation of intercepted water
            %             if (.not. use_vegetation.mlcanopyinst.h2ocan)
            %                 vegetation.mlcanopyinst.h2ocan(p,ic) = 0;
            %             end
            
%             vegetation.mlcanopyinst.h2ocan = h2ocan;  
        end  
    end
end
