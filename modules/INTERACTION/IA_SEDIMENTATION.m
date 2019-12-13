classdef IA_SEDIMENTATION

    properties
        IA_PARENT
        IA_CHILD
        STATUS
    end
    
    methods
        function ia_sedimentation = get_boundary_condition_u(ia_sedimentation, forcing)
            
           %get surfaceState
           ia_sedimentation.IA_CHILD.TEMP.surfaceState = forcing.TEMP.surfaceState;
           
           %if sedimentation is running and some sediment is there already,
           %get boundary conditions
           %if we run this with layerThick == 0, we get infinity flux.
           if ia_sedimentation.STATUS == 1 && ia_sedimentation.IA_CHILD.STATVAR.layerThick(1) > 0 
               child = ia_sedimentation.IA_CHILD;
               parent = ia_sedimentation.IA_PARENT;
               %call the native function for the sediment class 
               %this yields the same boundary conditions as the parent
               child = get_boundary_condition_u(child, forcing); 
               
               %calculate fluxes to parent and assign new boundary
               %condition (this is copied from IA_HEAT_SALT
               heatflux = (child.STATVAR.T(end) - parent.STATVAR.T(1)) .* child.STATVAR.thermCond(end) .* parent.STATVAR.thermCond(1) ./...
                    (child.STATVAR.thermCond(end).* parent.STATVAR.layerThick(1)./2 + parent.STATVAR.thermCond(1).* child.STATVAR.layerThick(end)./2 );
               child.TEMP.heatFlux_lb = -heatflux;
               parent.TEMP.heatFlux_ub = -heatflux;
               parent.TEMP.T_ub = child.STATVAR.T(end);
                           
               saltflux = (child.STATVAR.saltConc(end) - parent.STATVAR.saltConc(1)) .* child.STATVAR.saltDiff(end) .* parent.STATVAR.saltDiff(1) ./...
                    (child.STATVAR.saltDiff(end).* parent.STATVAR.layerThick(1)./2 +  parent.STATVAR.saltDiff(1).* child.STATVAR.layerThick(end)./2 );
               child.TEMP.saltFlux_lb = -saltflux;
               parent.TEMP.saltFlux_ub = -saltflux;
           end 
        end
        
        function ia_sedimentation = get_derivative_temperature_salt(ia_sedimentation)
            if ia_sedimentation.STATUS == 1 
                child = ia_sedimentation.IA_CHILD;
                child = get_derivatives_prognostic(child);               
            end
        end
        
        function ia_sedimentation = advance_prognostic(ia_sedimentation, timestep)
            if ia_sedimentation.STATUS == 1 
                child = ia_sedimentation.IA_CHILD;
                
                %% what needs to happen here?
                % advance in time
                % as a ghost cell
                % the same as for the parent?
                child = advance_prognostic(child, timestep);
            end
            if ia_sedimentation.IA_CHILD.TEMP.surfaceState == 1 || ia_sedimentation.IA_CHILD.TEMP.surfaceState == 0
                %calculate sedimentation
                %timestep is in days
                %sedRate is in m/year
                %accumulate sediment here, it gets set to zero when it is
                %sedimented in compute_diagnostic
                child = ia_sedimentation.IA_CHILD;
                child.TEMP.sedimentation = child.TEMP.sedimentation + child.PARA.sedRate * timestep/365.25;
            end
        end
        
        function ia_sedimentation = compute_diagnostic(ia_sedimentation)
            if ia_sedimentation.STATUS == 1 || ia_sedimentation.IA_CHILD.TEMP.surfaceState == 1 || ia_sedimentation.IA_CHILD.TEMP.surfaceState == 0
                
                %% what needs to happen here?
                % only update of thermal properties?
                % compaction?
                minLayerThick = 0.1;
                maxLayerThick = 0.25;
                
                child = ia_sedimentation.IA_CHILD;
                newLayerThick = child.TEMP.sedimentation;
                
                %if accumulated sediment is less than 10cm
                if newLayerThick < minLayerThick && child.STATVAR.layerThick(1) < minLayerThick 
                    %if upper layer of child sediment is less than
                    %10cm, and already accumulated sediment is also less 
                    %than 10cm, let sediment accumulate some more
                    
                    %this is the case for a new child
                    
                    %if we allow layers with less than 10cm, the
                    %timesteps are too small

                    %do nothing
                    
                elseif child.STATVAR.layerThick(1) < maxLayerThick
                    %if upper layer of sediment is more than 10cm but less
                    %than 50cm, add to upper layer
                    %accumulated sediment can be very small
                    
                    %if this layer was empty before, update temperature etc
                    if child.STATVAR.layerThick(1) == 0
                        parent = ia_sedimentation.IA_PARENT;
                        child.STATVAR.T = parent.STATVAR.T(1);
                        child.STATVAR.saltConc = parent.STATVAR.saltConc(1);
                        child.STATVAR.c_eff = parent.STATVAR.c_eff(1);
                        child.STATVAR.thermCond = parent.STATVAR.thermCond(1);
                        child.STATVAR.Tmelt = parent.STATVAR.Tmelt(1);
                        child.STATVAR.liqWater = parent.STATVAR.liqWater(1);
                    end
                                             
                        
                    child.STATVAR.layerThick(1) = child.STATVAR.layerThick(1)+newLayerThick;
                    
                    
                    %should we do something about conversion of energy?
                        
                    child.TEMP.sedimentation = 0; %reset sedimentation
                    
                elseif child.STATVAR.layerThick(1) >= maxLayerThick && newLayerThick > minLayerThick
                    %if upper layer is more than 50cm and accumulated
                    %sediment is more than 10cm, make new layer
                    child.STATVAR.layerThick = [newLayerThick; child.STATVAR.layerThick];

                    %expand state variables to new layer
                    child.STATVAR.T = [child.STATVAR.T(1); child.STATVAR.T];
                    child.STATVAR.saltConc = [child.STATVAR.saltConc(1); child.STATVAR.saltConc];
                    child.STATVAR.c_eff = [child.STATVAR.c_eff(1); child.STATVAR.c_eff];
                    child.STATVAR.thermCond = [child.STATVAR.thermCond(1); child.STATVAR.thermCond];
                    child.STATVAR.Tmelt = [child.STATVAR.Tmelt(1); child.STATVAR.Tmelt];
                    
                    %here maybe some compaction?
                    child.STATVAR.water = [child.STATVAR.water(1); child.STATVAR.water];
                    child.STATVAR.liqWater = [child.STATVAR.liqWater(1); child.STATVAR.liqWater];
                    child.STATVAR.porosity = [child.STATVAR.porosity(1); child.STATVAR.porosity];
                    child.STATVAR.mineral = [child.STATVAR.mineral(1); child.STATVAR.mineral];
                    child.STATVAR.organic = [child.STATVAR.organic(1); child.STATVAR.organic];

                    child.STATVAR.soilType = [child.STATVAR.soilType(1); child.STATVAR.soilType]; %is this a skalar or a vector?

                    child.TEMP.sedimentation = 0; %reset sedimentation
                else
                    %i.e. if upper layer is more than 50cm, but acumulated
                    %sediment is less than 10cm, let accumulate some more
                    
                    %nothing
                end
                    
                %update the thermal properties
                child = getThermalProps_wSalt(child);

                %update upper position
                child.STATVAR.upperPos = child.STATVAR.lowerPos + sum(child.STATVAR.layerThick);

                %set status to 1 (may already be so)
                ia_sedimentation.STATUS = 1;
            end
        end
        
       
        function ia_sedimentation = check_trigger(ia_sedimentation)
            % check if child is big enough to be own class
            Height = abs(ia_sedimentation.IA_CHILD.STATVAR.upperPos - ia_sedimentation.IA_CHILD.STATVAR.lowerPos);
           if ia_sedimentation.STATUS == 1 
               
               if Height > 2  || ia_sedimentation.IA_CHILD.STATVAR.soilType(1) ~= ia_sedimentation.IA_CHILD.TEMP.surfaceState
                   %make child the new top module
                   ia_sedimentation.IA_CHILD.PREVIOUS = ia_sedimentation.IA_PARENT.PREVIOUS;
                   ia_sedimentation.IA_CHILD.PREVIOUS.NEXT = ia_sedimentation.IA_CHILD;

                   ia_sedimentation.IA_CHILD.NEXT = ia_sedimentation.IA_PARENT;
                   ia_sedimentation.IA_PARENT.PREVIOUS = ia_sedimentation.IA_CHILD;

                   %add interaction
                   temp_store = ia_sedimentation.IA_CHILD; %the snow class

                   ia_class = get_IA_class(class(temp_store), class(temp_store.NEXT)); 
                   temp_store.IA_NEXT = ia_class;
                   temp_store.IA_NEXT.PREVIOUS = temp_store;
                   temp_store.IA_NEXT.NEXT = temp_store.NEXT;
                   temp_store.NEXT.IA_PREVIOUS = temp_store.IA_NEXT;

                   % destroy the class
                   ia_sedimentation = [];  
               end
           end
        end
        
    end
    
end

