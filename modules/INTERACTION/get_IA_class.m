function ia_class = get_IA_class(class1, class2) %change this later, class1 and 2 are strings


%ia_class = IA_HEAT_SALT;

%replace by matrix, this is stupid

if strcmp((class1), 'GROUND_freeW_seb') && strcmp((class2), 'GROUND_freeW_seb')
    ia_class = IA_HEAT();
elseif strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'GROUND_freeW_seb_snow')
    ia_class = IA_HEAT();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'GROUND_freeW_seb')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'GROUND_freeW_seb'))
    ia_class = IA_HEAT();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'SNOW_seb_simple')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'SNOW_seb_simple'))
    ia_class = IA_HEAT();
elseif (strcmp((class1), 'GROUND_freeW_seb') && strcmp((class2), 'SNOW_seb_simple')) || (strcmp((class2), 'GROUND_freeW_seb') && strcmp((class1), 'SNOW_seb_simple'))
    ia_class = IA_HEAT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class2), 'SNOW_seb_simple') || strcmp((class2), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class1), 'SNOW_seb_simple')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class2), 'GROUND_freeW_seb') || strcmp((class2), 'GROUND_fcSimple_salt_seb_snow') && strcmp((class1), 'GROUND_freeW_seb')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb') && strcmp((class2), 'GROUND_freeW_seb') || strcmp((class2), 'GROUND_fcSimple_salt_seb') && strcmp((class1), 'GROUND_freeW_seb')
    ia_class = IA_HEAT_SALT();
elseif strcmp((class1), 'GROUND_fcSimple_salt_seb') && strcmp((class2), 'GROUND_freeW_seb_snow') || strcmp((class2), 'GROUND_fcSimple_salt_seb') && strcmp((class1), 'GROUND_freeW_seb_snow')
    ia_class = IA_HEAT_SALT();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'GROUND_vegetation')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'GROUND_vegetation'))
    ia_class = IA_HEAT_VEGETATION();
elseif (strcmp((class1), 'GROUND_vegetation') && strcmp((class2), 'GROUND_freeW_seb')) || (strcmp((class2), 'GROUND_vegetation') && strcmp((class1), 'GROUND_freeW_seb'))
    ia_class = IA_HEAT_VEGETATION();
elseif (strcmp((class1), 'GROUND_freeW_seb_snow') && strcmp((class2), 'GROUND_vegetation_snow')) || (strcmp((class2), 'GROUND_freeW_seb_snow') && strcmp((class1), 'GROUND_vegetation_snow'))
    ia_class = IA_HEAT_VEGETATION();
    
elseif (strcmp((class1), 'GROUND_vegetation_snow') && strcmp((class2), 'GROUND_freeW_seb')) || (strcmp((class2), 'GROUND_vegetation_snow') && strcmp((class1), 'GROUND_freeW_seb'))  %%%% Used at the moment!!
    ia_class = IA_HEAT_VEGETATION();
elseif (strcmp((class1), 'GROUND_vegetation_snow') && strcmp((class2), 'SNOW_seb_simple_vegetation')) || (strcmp((class2), 'GROUND_vegetation_snow') && strcmp((class1), 'SNOW_seb_simple_vegetation'))  %%%% Used at the moment!!
    ia_class = IA_SNOW_GROUND_VEGETATION();
    
elseif (strcmp((class1), 'GROUND_vegetation_snow') && strcmp((class2), 'SNOW_seb_simple_vegetation_below')) || (strcmp((class2), 'GROUND_vegetation_snow') && strcmp((class1), 'SNOW_seb_simple_vegetation_below'))
    ia_class = IA_SNOW_GROUND_VEGETATION();
else
    ia_class=0;
    disp('combination of classes not supported')
end
end