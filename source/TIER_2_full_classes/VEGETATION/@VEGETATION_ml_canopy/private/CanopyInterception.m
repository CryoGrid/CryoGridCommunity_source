function vegetation = CanopyInterception(vegetation) %, p, ic, il

% %DESCRIPTION:
% Interception and throughfall

% p = vegetation.mlcanopyinst.p;

% %ARGUMENTS:
num_exposedvegp = vegetation.canopy.num_exposedvegp;       % Number of non-snow-covered veg points in CLM patch filter
filter_exposedvegp = vegetation.canopy.filter_exposedvegp;  % CLM patch filter for non-snow-covered vegetation
dtime_sub = vegetation.params.dtime_sub; % Model time step, for now, 1
% use_h2ocan = vegetation.physcon.use_h2ocan; % Canopy layer intercepted water (kg H2O/m2), for now, 1
%--------------------------------------------------------------------
% 1. CanopyWater
% 1.1 CanopyInterception
%Input
% lai(p)  = vegetation.canopy.lai(p);                                            % Leaf area index of canopy (m2/m2)
% sai  = vegetation.mlcanopyinst.sai;                                                  % Stem area index of canopy (m2/m2)
% ncan  = vegetation.mlcanopyinst.ncan;                                                % Number of aboveground layers
% dpai  = vegetation.canopy.dpai;                                                % Layer plant area index (m2/m2)
% dlai   = vegetation.canopy.dlai;                                               % Layer leaf area index (m2/m2)
% qflx_rain = vegetation.mlcanopyinst.qflx_rain; % forcing.DATA.rainfall               % Rainfall (mm H2O/s = kg H2O/m2/s)
% qflx_snow  = vegetation.mlcanopyinst.qflx_snow; % forcing.DATA.rainfall              % Snowfall (mm H2O/s = kg H2O/m2/s)

% %Input/Output
% h2ocan = vegetation.mlcanopyinst.h2ocan;                                             % Canopy layer intercepted water (kg H2O/m2)

%Output
% fwet = vegetation.mlcanopyinst.fwet;                                                 % Fraction of plant area index that is wet
% fdry = vegetation.mlcanopyinst.fdry;                                                 % Fraction of plant area index that is green and dry
% qflx_prec_intr = vegetation.mlcanopyinst.qflx_prec_intr;                             %Intercepted precipitation (kg H2O/m2/s)

%--------------------------------------------------------------------

% Time step
dtime = dtime_sub;

% Maximum allowed interception (kg H2O/m2 leaf)
dewmx = 0.1;

for f = 1:num_exposedvegp
    p = filter_exposedvegp(f);
    
    % Fraction of precipitation that is rain and snow
    if ((vegetation.mlcanopyinst.qflx_snow(p) + vegetation.mlcanopyinst.qflx_rain(p)) > 0)   %>0
        fracrain = vegetation.mlcanopyinst.qflx_rain(p) ./ (vegetation.mlcanopyinst.qflx_snow(p) + vegetation.mlcanopyinst.qflx_rain(p));
        fracsnow = vegetation.mlcanopyinst.qflx_snow(p) ./ (vegetation.mlcanopyinst.qflx_snow(p) + vegetation.mlcanopyinst.qflx_rain(p));
    else
        fracrain = 0;
        fracsnow = 0;
    end
    
    % Fraction of precipitation that is intercepted
    % fpi = 0.25_r8 * (1._r8 - exp(-0.5_r8*(lai(p) + sai(p))))   % CLM4.5
    interception_fraction = 0.5;
    fpi = interception_fraction .* tanh(vegetation.canopy.lai(p) + vegetation.mlcanopyinst.sai(p));        % CLM5
    
    % Direct throughfall
    qflx_through_rain = vegetation.mlcanopyinst.qflx_rain(p) .* (1. - fpi);
    qflx_through_snow = vegetation.mlcanopyinst.qflx_snow(p) .* (1. - fpi);
    
    % Intercepted precipitation
    vegetation.mlcanopyinst.qflx_prec_intr(p) = (vegetation.mlcanopyinst.qflx_snow(p) + vegetation.mlcanopyinst.qflx_rain(p)) .* fpi;
    
    % Find number of layers with lai+sai
    
    n = 0;
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
        if (vegetation.canopy.dpai(p,ic) > 0) %>0
            n = n + 1;
        end
    end
    
    % Loop through layers for water balance calculation
    qflx_candrip = 0;
    for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
        
        if (vegetation.canopy.dpai(p,ic) > 0)  % leaf layer
            
            % Maximum external water held in layer
            h2ocanmx = dewmx .* vegetation.canopy.dpai(p,ic);
            
            %SEBAS changed!!!!!!
            precip_fraction_per_layer = vegetation.canopy.dpai(p,ic) ./ vegetation.canopy.pai;
            %SEBAS:THIS PRODUCES by lai too high throughfall if the canopy
            %is completly wet - precip is applied to all layers  
            
            
            %             Water storage of intercepted precipitation. Intercepted water
            %             is applied equally to all layers.
            % % %             if use_h2ocan == 1
            
            %SEBAS changed!!!!!!
            %vegetation.mlcanopyinst.h2ocan(p,ic) = vegetation.mlcanopyinst.h2ocan(p,ic) + vegetation.mlcanopyinst.qflx_prec_intr(p) .* dtime; %dtime/float(n);
            vegetation.mlcanopyinst.h2ocan(p,ic) = vegetation.mlcanopyinst.h2ocan(p,ic) + precip_fraction_per_layer .* vegetation.mlcanopyinst.qflx_prec_intr(p) .* dtime; %dtime/float(n);
            %end SEBAS changed!!!!!!
            
            % % %             else
            % % %                 vegetation.mlcanopyinst.h2ocan(p,ic) = 0;
            % % %             end
            
            % Excess water that exceeds the maximum capacity. If xrun > 0
            % then h2ocan is set to h2ocanmx and excess water is added to
            % throughfall.

            xrun = (vegetation.mlcanopyinst.h2ocan(p,ic) - h2ocanmx) / dtime;
            if (xrun > 0.)
               
                qflx_candrip = qflx_candrip + xrun;
                vegetation.mlcanopyinst.h2ocan(p,ic) = h2ocanmx;
            end
            
            % Wetted fraction of canopy
            maximum_leaf_wetted_fraction = 0.05;
            vegetation.mlcanopyinst.fwet(p,ic) = max((vegetation.mlcanopyinst.h2ocan(p,ic)/h2ocanmx),0.)^0.67;
            vegetation.mlcanopyinst.fwet(p,ic) = min(vegetation.mlcanopyinst.fwet(p,ic), maximum_leaf_wetted_fraction);
            
            % Fraction of canopy that is green and dry
            vegetation.mlcanopyinst.fdry(p,ic) = (1 - vegetation.mlcanopyinst.fwet(p,ic)) .* (vegetation.canopy.dlai(p,ic) ./ vegetation.canopy.dpai(p,ic));
            
        else % non-leaf layer
            vegetation.mlcanopyinst.h2ocan(p,ic) = 0.;
            vegetation.mlcanopyinst.fwet(p,ic) = 0.;
            vegetation.mlcanopyinst.fdry(p,ic) = 0.;
            
        end
        
    end
    
    % Total throughfall onto ground
    vegetation.mlcanopyinst.qflx_prec_grnd_rain = qflx_through_rain + qflx_candrip .* fracrain; % mm h2o/s
    vegetation.mlcanopyinst.qflx_prec_grnd_snow = qflx_through_snow + qflx_candrip .* fracsnow;

    %     disp([vegetation.mlcanopyinst.qflx_rain vegetation.mlcanopyinst.qflx_prec_grnd_rain])
    

    
    % Output
    % vegetation.mlcanopyinst.h2ocan = h2ocan;                                  %=> mlcanopy_inst%h2ocan            % Canopy layer intercepted water (kg H2O/m2)
    % vegetation.mlcanopyinst.fwet = fwet;                                      %=> mlcanopy_inst%fwet              % Fraction of plant area index that is wet
    % vegetation.mlcanopyinst.fdry = fdry;                                      %=> mlcanopy_inst%fdry              % Fraction of plant area index that is green and dry
    % vegetation.mlcanopyinst.qflx_prec_intr = qflx_prec_intr;                  %=> mlcanopy_inst%qflx_prec_intr    % Intercepted precipitation (kg H2O/m2/s)
    
    
end
end

