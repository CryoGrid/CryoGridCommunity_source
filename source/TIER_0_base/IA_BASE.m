%========================================================================
% CryoGrid base class for all INTERACTION (IA) classes 
% S. Westermann, Oct 2020
%========================================================================

classdef IA_BASE < matlab.mixin.Copyable 

     properties
        PREVIOUS
        NEXT
     end
     
     methods
         function finalize_init(IA_BASE, tile) 
            %do nothing
         end
     end
end