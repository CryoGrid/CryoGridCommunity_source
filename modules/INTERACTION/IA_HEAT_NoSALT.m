%> IA_HEAT_SALT heat interaction with zero-salt-flux between two classes
%> transition from salty to non-salty class
%> heat is transfered between two classes via a heatflux proportional to
%> the temperature difference at the border
%> salt flux is set to zero
classdef IA_HEAT_NOSALT < IA_HEAT 
    %transition from salty to non-salty class, zero flux boundary
    %condition imposed at both boundaries, this should be changd depnding on context of use

    

    methods
        
        function get_boundary_condition_m(ia_heat_salt)
            % calculate heat flux with the super class method
            get_boundary_condition_m@IA_HEAT(ia_heat_salt)
            
            % get the two classes involved
            stratigraphy1 = ia_heat_salt.PREVIOUS;
            stratigraphy2 = ia_heat_salt.NEXT;
            
            % set the salt flux to zero
            stratigraphy1.TEMP.saltConcFlux_lb = 0;
            stratigraphy2.TEMP.saltConcFlux_ub = 0; 
        end
    end
end
   
   
