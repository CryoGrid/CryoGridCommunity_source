%> IA_HEAT_SALT_FLUX heat and salt interaction between two classes
%> heat is transfered between two classes via a heatflux proportional to
%> the temperature difference at the border
%> salt is transferred between two classes via a saltflux proportional to
%> the salt concentration difference at the border
classdef IA_HEAT_SALT_FLUX < IA_HEAT
    
    methods
        
        function get_boundary_condition_m(ia_heat_salt_flux)
            % calculate heat flux with the super class method
            get_boundary_condition_m@IA_HEAT(ia_heat_salt_flux)
            
            % get the two classes involved
            stratigraphy1 = ia_heat_salt_flux.PREVIOUS;
            stratigraphy2 = ia_heat_salt_flux.NEXT;
            
            % claculate salt flux based on the salt concentration
            % difference at the border
            saltflux = (stratigraphy1.STATVAR.saltConc(end) - stratigraphy2.STATVAR.saltConc(1)) .* stratigraphy1.STATVAR.saltDiff(end) .* stratigraphy2.STATVAR.saltDiff(1) ./...
                (stratigraphy1.STATVAR.saltDiff(end).* stratigraphy2.STATVAR.layerThick(1)./2 +  stratigraphy2.STATVAR.saltDiff(1).* stratigraphy1.STATVAR.layerThick(end)./2 );
            
            % assign the flux to the lower/upper boundary of the
            % upper/lower class, respectively
            stratigraphy1.TEMP.saltFlux_lb = -saltflux;
            stratigraphy2.TEMP.saltFlux_ub = -saltflux;
            
        end
    end
end
   