%========================================================================
% CryoGrid get_IA_class function
% identifies and initializes the correct interaction class between pairs of
% GROUND classes, returns 0 if classes are not compatible
% INPUT: GROUND class names as strings
% NOTE: be very careful when updating the interaction matrix "ia"! Keep the
% last working version in the code (commented out)
% S. Westermann, October 2020
% modified to make code additive
% S. Westermann, December 2021
%========================================================================


function ia_class = get_IA_class(above_class, below_class)


ia_class = 0;

%modfications of class strings

%FLIP_FLOP classses have same interactions as original classes
if size(above_class,2) > 10 && strcmp(above_class(end-8:end), 'FLIP_FLOP')
    above_class = above_class(1:end-10);
end

if size(below_class,2) > 10 && strcmp(below_class(end-8:end), 'FLIP_FLOP')
    below_class = below_class(1:end-10);
end


if strcmp(above_class, 'Top') %allows treating Top like a "normal" class, as IA class is never called between Top and the following class 
    ia_class = IA_BASE();
    
elseif strcmp(above_class, 'GROUND_freeW_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11();
    end
    if strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER01();
    end
    if strcmp(below_class, 'GROUND_fcSimple_salt_seb') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
elseif strcmp(above_class, 'GROUND_freeW_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER01();
        
    elseif strcmp(below_class, 'GROUND_fcSimple_salt_seb') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
elseif strcmp(above_class, 'GROUND_freeW_bucketW_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11();
    end
    
elseif strcmp(above_class, 'GROUND_freeW_bucketW_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb')
        ia_class = IA_HEAT11_WATER11();
    end
    
elseif strcmp(above_class, 'GROUND_freeW_bucketW_convection_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11();
    end
    
elseif strcmp(above_class, 'GROUND_freeW_bucketW_convection_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb')
        ia_class = IA_HEAT11_WATER11();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER01();
        
    elseif strcmp(below_class, 'GROUND_fcSimple_salt_seb') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER01();
        
    elseif strcmp(below_class, 'GROUND_fcSimple_salt_seb') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_bucketW_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_bucketW_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb')
        ia_class = IA_HEAT11_WATER11();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_bucketW_Xice_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(above_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_RichardsEqW_seb') || strcmp(above_class, 'GROUND_freezeC_RichardsEqW_seb_pressure')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_RichardsEqW_seb_snow')  || strcmp(above_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_RichardsEqW_Xice_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
    end
    
elseif strcmp(above_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_WATER10();
    end
    
elseif strcmp(above_class, 'GROUND_fcSimple_salt_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_SALT10();
    end
    
elseif strcmp(above_class, 'GROUND_fcSimple_salt_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11_SALT10();
    end
    
elseif strcmp(above_class, 'LAKE_simple_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'LAKE_simple_unfrozen_seb')
        ia_class = IA_LAKE_simple_frozen_unfrozen();
    end
    
elseif strcmp(above_class, 'LAKE_simple_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'LAKE_simple_unfrozen_seb')
        ia_class = IA_LAKE_simple_frozen_unfrozen();
    end
    
elseif strcmp(above_class, 'LAKE_simple_unfrozen_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11_LAKE();
        
    elseif strcmp(below_class, 'LAKE_simple_seb') || strcmp(below_class, 'LAKE_simple_seb_snow')
        ia_class = IA_LAKE_simple_frozen_unfrozen();
    end
    
elseif strcmp(above_class, 'LAKE_simple_bucketW_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11_WATER10_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER11_LAKE_XICE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_LAKE_XICE();
        
    elseif strcmp(below_class, 'LAKE_simple_bucketW_seb') || strcmp(below_class, 'LAKE_simple_bucketW_seb_snow') || strcmp(below_class, 'LAKE_simple_unfrozen_bucketW_seb')
        ia_class = IA_LAKE_simple_frozen_unfrozen();
    end
    
elseif strcmp(above_class, 'LAKE_simple_bucketW_seb_snow')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11_WATER10_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER11_LAKE_XICE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_LAKE_XICE();
        
    elseif strcmp(below_class, 'LAKE_simple_bucketW_seb') || strcmp(below_class, 'LAKE_simple_bucketW_seb_snow') || strcmp(below_class, 'LAKE_simple_unfrozen_bucketW_seb')
        ia_class = IA_LAKE_simple_frozen_unfrozen();
    end
    
elseif strcmp(above_class, 'LAKE_simple_unfrozen_bucketW_seb')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11_WATER10_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER11_LAKE_XICE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_LAKE_XICE();
        
    elseif strcmp(below_class, 'LAKE_simple_bucketW_seb') || strcmp(below_class, 'LAKE_simple_bucketW_seb_snow') || strcmp(below_class, 'LAKE_simple_unfrozen_bucketW_seb')
        ia_class = IA_LAKE_simple_frozen_unfrozen();
    end
    
    
     %TTOP as upper boundary condition
elseif strcmp(above_class, 'GROUND_TTOP_simple2')
    if strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')  || strcmp(below_class, 'GROUND_freezeC_seb') ...
            || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') ...
            || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb')
        ia_class = IA_HEAT_TTOP();
    end
        
    %SNOW classes
elseif strcmp(above_class, 'SNOW_simple_seb')
    if strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb_snow') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow') || strcmp(below_class, 'LAKE_simple_seb_snow') || strcmp(below_class, 'GLACIER_freeW_seb_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow') || strcmp(below_class, 'LAKE_simple_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER01();
    end
    
elseif strcmp(above_class, 'SNOW_simple_bucketW_seb')
    if strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb_snow') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow') || strcmp(below_class, 'LAKE_simple_seb_snow')  || strcmp(below_class, 'GLACIER_freeW_seb_snow')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER11_SNOW_XICE();
        
    elseif strcmp(below_class, 'LAKE_simple_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_SNOW_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW_XICE();
    end
    
elseif strcmp(above_class, 'SNOW_crocus_bucketW_seb')
    if strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb_snow') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow') || strcmp(below_class, 'LAKE_simple_seb_snow') || strcmp(below_class, 'GLACIER_freeW_seb_snow')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER11_SNOW_XICE();
        
    elseif strcmp(below_class, 'LAKE_simple_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_SNOW_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW_XICE();
    end
    
elseif strcmp(above_class, 'SNOW_crocus_bucketW_seb_vegetation')
    if strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb_snow') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow') || strcmp(below_class, 'LAKE_simple_seb_snow')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER11_SNOW_XICE();
        
    elseif strcmp(below_class, 'LAKE_simple_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_SNOW_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW_XICE();
    
    elseif strcmp(below_class, 'VEGETATION_CLM5_seb_snow')
        ia_class = IA_SEB_vegetation_CLM5_SNOW();
    end
    
elseif strcmp(above_class, 'SNOW_crocus2_bucketW_seb')
    if strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb_snow') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow') || strcmp(below_class, 'LAKE_simple_seb_snow') || strcmp(below_class, 'GLACIER_freeW_seb_snow')
        ia_class = IA_HEAT11_WATER10();
        
    elseif strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER11_SNOW_XICE();
        
    elseif strcmp(below_class, 'LAKE_simple_bucketW_seb_snow')
        ia_class = IA_HEAT11_WATER11_SNOW_LAKE();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_pressure_snow') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_vegetation_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow')
        ia_class = IA_HEAT11_WATER11_RichardsEq_SNOW_XICE();
    
    elseif strcmp(below_class, 'VEGETATION_CLM5_seb_snow')
        ia_class = IA_SEB_vegetation_CLM5_SNOW();
    end
    
    %-----------store_flip_flop -----
    
elseif strcmp(above_class, 'GROUND_store_flip_flop_singleClass') || strcmp(above_class, 'GROUND_store_flip_flop_singleClass_BGC')
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER01();
        
    elseif strcmp(below_class, 'GROUND_fcSimple_salt_seb') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
    %------------GLACIER----
    
elseif strcmp(above_class, 'GLACIER_freeW_seb')  %modify if more interactions for water flow are possible!!
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freezeC_seb')
        ia_class = IA_HEAT11();
    end
    if strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER01();
    end
    if strcmp(below_class, 'GROUND_fcSimple_salt_seb') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
elseif strcmp(above_class, 'GLACIER_freeW_seb_snow') %modify if more interactions for water flow are possible!!
    if strcmp(below_class, 'GROUND_freeW_seb') || strcmp(below_class, 'GROUND_freeW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_seb') || strcmp(below_class, 'GROUND_freezeC_seb_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_freeW_bucketW_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb') || strcmp(below_class, 'GROUND_freeW_bucketW_convection_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow') || strcmp(below_class, 'GROUND_freezeC_bucketW_Xice_seb_snow_BGC')
        ia_class = IA_HEAT11_WATER01();
        
    elseif strcmp(below_class, 'GROUND_fcSimple_salt_seb') || strcmp(below_class, 'GROUND_fcSimple_salt_seb_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
        % ---------- VEGETATION -----------
    
elseif strcmp(above_class, 'VEGETATION_CLM5_seb')
    if strcmp(below_class, 'SNOW_crocus_bucketW_seb') || strcmp(below_class, 'SNOW_simple_bucketW_seb') || strcmp(below_class, 'SNOW_crocus2_bucketW_seb')
        ia_class = IA_SEB_vegetation_CLM5_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb') || strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb')
        ia_class = IA_SEB_vegetation_CLM5();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow')
        ia_class = IA_SEB_vegetation_CLM5_GROUND_snow();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow') 
        ia_class = IA_SEB_vegetation_CLM5_GROUND_Xice_snow();
    end
    
       % ---------- VEGETATION with SNOW in canopy-----------
    
elseif strcmp(above_class, 'VEGETATION_CLM5_seb_snow')
    if strcmp(below_class, 'SNOW_crocus_bucketW_seb') || strcmp(below_class, 'SNOW_crocus2_bucketW_seb')
        ia_class = IA_SEB_vegetation_CLM5_SNOW();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_seb_snow')
        ia_class = IA_SEB_vegetation_CLM5_GROUND_snow();
        
    elseif strcmp(below_class, 'GROUND_freezeC_RichardsEqW_Xice_seb_snow') 
        ia_class = IA_SEB_vegetation_CLM5_snow_GROUND_Xice_snow();
    end
    
    %--------MULTI-TILE----
    
elseif strcmp(above_class, 'GROUND_multi_tile2') || strcmp(above_class, 'GROUND_multi_tile2_Cas')
    ia_class = IA_MULTI_TILE10();
    
    
    
    %----Dirichlet BC classes-----
    
elseif strcmp(above_class, 'GROUND_freeW_ubT') || strcmp(above_class, 'GROUND_freeW_ubT_snow')
    if strcmp(below_class, 'GROUND_freezeC_ubT') || strcmp(below_class, 'GROUND_freezeC_ubT_snow') || strcmp(below_class, 'GROUND_freeW_ubT') || strcmp(below_class, 'GROUND_freeW_ubT_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_fcSimple_salt_ubT') || strcmp(below_class, 'GROUND_fcSimple_salt_ubT_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
    
elseif strcmp(above_class, 'GROUND_freezeC_ubT') || strcmp(above_class, 'GROUND_freezeC_ubT_snow')
    if strcmp(below_class, 'GROUND_freezeC_ubT') || strcmp(below_class, 'GROUND_freezeC_ubT_snow') || strcmp(below_class, 'GROUND_freeW_ubT') || strcmp(below_class, 'GROUND_freeW_ubT_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_fcSimple_salt_ubT') || strcmp(below_class, 'GROUND_fcSimple_salt_ubT_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
    
elseif strcmp(above_class, 'GROUND_fcSimple_salt_ubT') || strcmp(above_class, 'GROUND_fcSimple_salt_ubT_snow')
    if strcmp(below_class, 'GROUND_freezeC_ubT') || strcmp(below_class, 'GROUND_freezeC_ubT_snow') || strcmp(below_class, 'GROUND_freeW_ubT') || strcmp(below_class, 'GROUND_freeW_ubT_snow')
        ia_class = IA_HEAT11_SALT10();
    end
    
    
elseif strcmp(above_class, 'SNOW_simple_ubT')
    if strcmp(below_class, 'GROUND_freeW_ubT_snow') || strcmp(below_class, 'GROUND_freezeC_ubT_snow')
        ia_class = IA_HEAT11();
        
    elseif strcmp(below_class, 'GROUND_fcSimple_salt_ubT_snow')
        ia_class = IA_HEAT11_SALT01();
    end
    
elseif strcmp(above_class, 'GROUND_house')
    ia_class = IA_HEAT11();
    
elseif strcmp(above_class, 'GROUND_house_on_poles')
    ia_class = IA_thermRadiation11();
    
    
end

if ia_class == 0
    switch above_class
        case {'GROUND_fcSimple_salt_ubtf', 'GROUND_fcSimple_salt_ubtf_snow',...
                'GROUND_fcSimple_salt_ubT', 'GROUND_fcSimple_salt_ubT_snow',...
                'GROUND_fcSimple_salt_seb', 'GROUND_fcSimple_salt_seb_snow'}
            switch below_class
                case {'GROUND_fcSimple_salt_ubtf',...
                        'GROUND_fcSimple_salt_ubT',...
                        'GROUND_fcSimple_salt_seb'}
                    error('Class interaction not currently supported!')
                case {'GROUND_freeW_seb',...
                        'GROUND_freeW_ubT',...
                        'GROUND_freeW_ubtf',...
                        'GROUND_freezeC_seb',...
                        'GROUND_freezeC_ubT'}
                    ia_class = IA_HEAT11_SALT10();
            end
        case {'GROUND_freeW_seb', 'GROUND_freeW_seb_snow',...
                'GROUND_freeW_ubT', 'GROUND_freeW_ubT_snow',...
                'GROUND_freeW_ubtf', 'GROUND_freeW_ubtf_snow',...
                'GROUND_freezeC_seb', 'GROUND_freezeC_seb_snow',...
                'GROUND_freezeC_ubT', 'GROUND_freezeC_ubT_snow'}
            switch below_class
                case {'GROUND_fcSimple_salt_ubtf',...
                        'GROUND_fcSimple_salt_ubT',...
                        'GROUND_fcSimple_salt_seb'}
                    ia_class = IA_HEAT11_SALT01();
            end
        case {'SNOW_simple_seb',...
                'SNOW_simple_ubtf_mf',...
                'SNOW_simple_ubT'}
            switch below_class
                case {'GROUND_fcSimple_salt_ubtf',...
                        'GROUND_fcSimple_salt_ubT',...
                        'GROUND_fcSimple_salt_seb'}
                    ia_class = IA_HEAT11_SALT01();
                case {'GROUND_freeW_ubtf_snow'}
                    ia_class = IA_HEAT11();
            end
    end
    
    
end
end


